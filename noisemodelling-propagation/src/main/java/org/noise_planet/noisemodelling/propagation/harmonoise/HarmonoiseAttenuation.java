/**
 * NoiseModelling is a library capable of producing noise maps. It can be freely used either for research and education, as well as by experts in a professional use.
 * <p>
 * NoiseModelling is distributed under GPL 3 license. You can read a copy of this License in the file LICENCE provided with this software.
 * <p>
 * Official webpage : http://noise-planet.org/noisemodelling.html
 * Contact: contact@noise-planet.org
 */

package org.noise_planet.noisemodelling.propagation.harmonoise;

import org.apache.commons.math3.complex.Complex;
import org.locationtech.jts.algorithm.Angle;
import org.locationtech.jts.geom.Coordinate;
import org.noise_planet.noisemodelling.propagation.SceneWithAttenuation;

import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

/**
 * Compute excess attenuation according to Harmonoise propagation model
 * Ref: Salomons, E., Van Maercke, D., Defrance, J.,&amp;De Roo, F. (2011). The Harmonoise sound propagation model.
 * Acta acustica united with acustica, 97(1), 62-74
 * @author Martin Glesser
 */

public class HarmonoiseAttenuation {
    SceneWithAttenuation scene; // Scene with attenuation data
    HarmonoiseAttenuationOutput attenuationOutput; // Output of the attenuation computation
    HarmonoiseGroundProfile groundProfile;
    private double[] waveNumber;

    // To be exposed to user (input of the model)
    double stdHs = 0; // standard deviation on the source height
    double stdHr = 0; // standard deviation on the receiver height
    double gammaT = 5 * Math.pow(10,-6); // typical value for moderate turbulence

    public HarmonoiseAttenuation(SceneWithAttenuation scene, HarmonoiseAttenuationOutput output) {
        this.scene = scene;
        this.attenuationOutput = output;
        this.groundProfile = new HarmonoiseGroundProfile(output.getCutProfile());
        waveNumber = scene.defaultCnossosParameters.getFrequenciesExact()
                .stream()
                .mapToDouble(f -> 2 * Math.PI * f / scene.defaultCnossosParameters.getCelerity())
                .toArray();
    }

    /**
     * Recursive calculation scheme for excess attenuation
     * Ref: section 2.2.3 from Salomons et al.
     */
    public void computeExcessAttenuation() {
        computeExcessAttenuation(0, groundProfile.getNVertices()-1);
    }

    /**
     * Recursive calculation scheme for excess attenuation
     * Ref: section 2.2.3 from Salomons et al.
     *
     * @param iStart index of the first ground profile vertex
     * @param iEnd index of the last ground profile vertex
     */
    private void computeExcessAttenuation(int iStart, int iEnd){
        Coordinate[] vertices = groundProfile.getVertices(iStart, iEnd);
        double maxDistance = 0;
        int indexMaxDistance = 0;
        for (int i = 1; i < vertices.length - 1; i++) {
            double crossProduct = (vertices[i].x - vertices[0].x) * (vertices[0].y - vertices[vertices.length-1].y)
                    + (vertices[i].y - vertices[0].y) * (vertices[vertices.length-1].x - vertices[0].x);
            if (crossProduct > 0) { // if the current point is above the line crossing the first and last points
                double distance = vertices[0].distance(vertices[i]) + vertices[i].distance(vertices[vertices.length-1])
                        - vertices[0].distance(vertices[vertices.length-1]); // Eq. 5
                if (distance > maxDistance){
                    maxDistance = distance;
                    indexMaxDistance = i;
                }
            }
        }
        if (indexMaxDistance == 0) { // No diffraction point above the line crossing the first and last points
            computeGroundAttenuation(iStart, iEnd);
        }else {
            computeDiffractionAttenuation(vertices[0], vertices[vertices.length-1],
                    vertices[indexMaxDistance]);
            computeExcessAttenuation(iStart, iStart + indexMaxDistance);
            computeExcessAttenuation(iStart + indexMaxDistance, iEnd);
        }
    }

    /**
     * Compute diffraction attenuation according to Harmonoise methodology and add it to AttenuationOutput
     * Ref: section 2.3 from Salomons et al.
     *
     * @param source source point ("real" or secondary at diffraction edge)
     * @param receiver receiver point ("real" or secondary at diffraction edge)
     * @param point diffraction point
     */
    private void computeDiffractionAttenuation(Coordinate source, Coordinate receiver, Coordinate point) {
        computeDiffractionAttenuation(source, receiver, point, new Complex[0]);
    }

    /**
     * Compute diffraction attenuation and diffracted sound pressure amplitude according to Harmonoise methodology
     * Ref: section 2.3 from Salomons et al.
     *
     * @param source source point ("real" or secondary at diffraction edge)
     * @param receiver receiver point ("real" or secondary at diffraction edge)
     * @param point diffraction point
     * @param diffractedPressure store the sound pressure if input size != 0
     */
    private void computeDiffractionAttenuation(Coordinate source, Coordinate receiver, Coordinate point, Complex[] diffractedPressure){
        double sourceAngle = - (Angle.angle(point, source) - Angle.PI_OVER_2);
        double receiverAngle = Angle.angle(point, receiver) + Angle.PI_OVER_2;
        double theta = sourceAngle + receiverAngle;

        double pathLengthDiff;
        double directPathLength;
        double distanceToSource = source.distance(point);
        double distanceToReceiver = receiver.distance(point);
        if (theta <= Math.PI) {
            directPathLength = Math.sqrt(Math.pow(distanceToSource, 2)
                    + Math.pow(distanceToReceiver, 2)
                    - 2 * distanceToSource * distanceToReceiver * Math.cos(theta)); // Eq. 10
            pathLengthDiff = - (distanceToSource + distanceToReceiver - directPathLength); // Eq. 9
        } else {
            directPathLength = distanceToSource + distanceToReceiver; // Eq. 12
            double epsilon = Math.sqrt(distanceToSource * distanceToReceiver) / (distanceToSource + distanceToReceiver)
                    * (theta - Math.PI); // Eq. 13
            pathLengthDiff = directPathLength * (1.0/2 * Math.pow(epsilon,2) + 1.0/3 * Math.pow(epsilon,3)); // Eq. 11
        }
        double[] fresnelNumber = scene.defaultCnossosParameters.getFrequenciesExact()
                .stream()
                .mapToDouble(f -> 2 * pathLengthDiff / (scene.defaultCnossosParameters.getCelerity() / f))
                .toArray(); // Eq.8
        double[] diffractionAttenuation = Arrays.stream(fresnelNumber)
                .map(HarmonoiseAttenuation::fresnelApproximation)
                .toArray();
        if (diffractedPressure.length == 0) {
            attenuationOutput.addDiffractionAttenuation(diffractionAttenuation);
        } else {
            Complex[] sphericalPressure = unitSphericalWavePressure(directPathLength);
            for (int i = 0; i < diffractionAttenuation.length; i++) {
                diffractedPressure[i] = sphericalPressure[i].multiply(10 * Math.log10(diffractionAttenuation[i] / 20));
            }
        }
    }

    /**
     * Return the sound pressure of unit-amplitude harmonic spherical wave
     * p(r) = exp(ikr)/r
     *
     * @param distance radial distance r from the source [m]
     * @return sound pressure p [Pa]
     */
    private Complex[] unitSphericalWavePressure(double distance){
        return (Complex[]) Arrays.stream(waveNumber)
                .mapToObj(k -> (new Complex(0,k * distance)).exp().divide(distance))
                .toArray();
    }

    /**
     * Approximation of the diffraction attenuation depending on the Fresnel number
     * Ref: Eq. 7 from Salomons et al.
     *
     * @param fresnelNumber Fresnel Number
     * @return diffraction attenuation
     */
    private static double fresnelApproximation(double fresnelNumber){
        if (fresnelNumber < -0.25) {
            return 0;
        } else if (fresnelNumber < 0) {
            return -6 + 12 * Math.sqrt(-fresnelNumber);
        } else if (fresnelNumber < 0.25) {
            return -6 - 12 * Math.sqrt(fresnelNumber);
        } else if (fresnelNumber < 1) {
            return -8 - 8 * Math.sqrt(fresnelNumber);
        } else {
            return -16 - 10 * Math.log10(fresnelNumber);
        }
    }

    private void computeGroundAttenuation(int iStart, int iEnd) {
        if (hasConvexSegment(groundProfile.getVertices(iStart, iEnd))){
            attenuationOutput.excessAttenuation += 0;
        }
        attenuationOutput.excessAttenuation += 0;
    }

    /**
     * Determine whether the ground profile contains convex segment or not
     * Ref: Figure 3 and "Geometry" subsection of section 2.4.1 from Salomons et al.
     *
     * @param vertices vertices of the ground profile
     * @return true if the profile contains at least one convex segment
     */
    private static boolean hasConvexSegment(Coordinate[] vertices){
        int iEnd = vertices.length - 1;
        boolean isConvex = false;
        // Loop on segments
        for (int i = 0; i < vertices.length - 2; i++) {
            double localSourceHeight = vertices[i+1].distance(vertices[0])
                    * Math.sin(Angle.angleBetweenOriented(vertices[0], vertices[i+1], vertices[i]));
            double localReceiverHeight = vertices[i].distance(vertices[iEnd])
                    * Math.sin(Angle.angleBetweenOriented(vertices[i+1], vertices[i], vertices[iEnd]));
            if (localSourceHeight < 0 || localReceiverHeight < 0){
                isConvex = true;
                break;
            }
        }
        return isConvex;
    }

    /**
     * Compute the spherical wave reflection coefficient
     * Ref: "Spherical-wave reflection coefficient" subsection of section 2.4.1 from Salomons et al.
     * Ref: K. Attenborough, K. M. Li, and K. Horoshenkov, Predicting Outdoor Sound. Taylor & Francis, 2006.
     * doi: 10.1201/9781482295023. (section 2.3)
     *
     * @param groundImpedance normalized ground impedance[]
     * @param angle reflection angle with respect to the normal on the segment
     * @param frequency frequency [Hz]
     * @param distance total distance between the image source and the receiver, through the segment
     * @param hm (hS + hR) / 2
     * @return spherical wave reflection coefficient
     */
    private Complex sphericalWaveReflectionCoefficient(Complex groundImpedance, double angle,
                                                             double frequency, double distance, double hm){
        Complex z =  groundImpedance.multiply(Math.cos(angle));
        Complex planeWaveReflectionCoefficient = z.subtract(1).divide(z.add(1));
        double kr = 2 * Math.PI * frequency / scene.defaultCnossosParameters.getCelerity() * distance;
        Complex admittance = groundImpedance.reciprocal();
        Complex numericalDistance = admittance.multiply(Math.sqrt(kr))
                .multiply(new Complex(0.5, 0.5));
        double hG = scene.defaultCnossosParameters.celerity / frequency / 32;
        double nG = 1 - 0.7 * Math.exp(-hm / hG);
        return boundaryLossFactor(numericalDistance).pow(nG)
                .multiply(new Complex(1).subtract(planeWaveReflectionCoefficient))
                .add(planeWaveReflectionCoefficient);
    }

    /**
     * Compute the boundary loss factor F(w)
     * Ref: K. Attenborough, K. M. Li, and K. Horoshenkov, Predicting Outdoor Sound. Taylor & Francis, 2006.
     * doi: 10.1201/9781482295023. (section 2.3)
     *
     * @param w numerical distance
     * @return boundary loss factor
     */
    private Complex boundaryLossFactor(Complex w){
        double x = w.getReal();
        double y = w.getImaginary();
        // Compute z = exp(-w^2) * erfc(-iw) for different value of the numerical distance
        Complex z;
        if ( Math.abs(x) > 3.9 || Math.abs(y) > 3){
            if ( Math.abs(x) > 6. || Math.abs(y) > 6.) {
                z = new Complex(0, 1).multiply(w)
                        .multiply(
                                w.pow(2).subtract(0.2752551).reciprocal().multiply(0.5124242)
                                .add(w.pow(2).subtract(2.724745).reciprocal().multiply(0.05176536))
                        );
            } else {
                z = new Complex(0, 1).multiply(w)
                        .multiply(
                                w.pow(2).subtract(0.1901635).reciprocal().multiply(0.461313500)
                                .add(w.pow(2).subtract(1.7844927).reciprocal().multiply(0.099992160))
                                .add(w.pow(2).subtract(5.5253437).reciprocal().multiply(0.002883894))
                        );
            }
        } else {
            double h = 0.8;
            double a1 = Math.cos(2 * x * y);
            double b1 = Math.sin(2 * x * y);
            double c1 = Math.exp(-2 * y * Math.PI / h) - Math.cos(2 * x * Math.PI / h);
            double d1 = Math.sin(2 * x * Math.PI / h);
            double cd = c1*c1 + d1*d1;
            double p = 1;
            double q = 1;
            if (cd != 0){
                double expArgument = -(Math.pow(x, 2) + 2 * y * Math.PI / h - Math.pow(y, 2));
                p = 2 * Math.exp(expArgument) * (a1 * c1 - b1 * d1) / cd;
                q = 2 * Math.exp(expArgument) * (a1 * d1 - b1 * c1) / cd;
            }
            double eh = Math.pow(10, -6);
            double h1 = 0;
            double h2 = 0;
            for (int n=1; n<=5; n++){
                double x1 = (y*y+x*x+n*n*h*h);
                double x2 = (y*y-x*x+n*n*h*h);
                double x3 = (y*y+x*x-n*n*h*h);
                double denominator = x2 * x2 + 4 * y * y * x * x;
                double numeratorFactor = Math.exp(-1 * (n * n) * h * h);
                h1 += numeratorFactor * x1 / denominator;
                h2 += numeratorFactor * x3 / denominator;
            }
            double denominator = Math.PI * (x * x + y * y);
            double k1 = h*y / denominator + 2*y*h*h1/Math.PI - y*eh/Math.PI;
            double k2 = h*x / denominator + 2*x*h*h2/Math.PI + x*eh/Math.PI;
            if (y < Math.PI/h)
            {
                k1 += p;
                k2 -= q;
            }
            if (y == Math.PI/h)
            {
                k1 += 0.5 * p;
                k2 -= 0.5 * q;
            }
            z = new Complex (k1, k2) ;
        }
        // return boundary loss factor
        return w.multiply(new Complex(0,1).multiply(Math.sqrt(Math.PI))).multiply(z).add(1);
    }

    /**
     * Compute the geometrical weighting factor corresponding to a segment
     * Ref: "Geometrical weighting factor" subsection of section 2.4.1 from Salomons et al.
     *
     * @param iSeg index of the first index of the segment
     * @param iStart first index of the ground profile or sub-profile
     * @param iEnd last index of the ground profile or sub-profile
     * @return geometrical weighting factor
     */
    private Complex[] geometricalWeightingFactor(int iSeg, int iStart, int iEnd){
        int nFreq = scene.defaultCnossosParameters.getFrequenciesExact().size();
        int lastVertex = groundProfile.getNVertices()-1;
        Coordinate realSource = groundProfile.getVertex(0);
        Coordinate realSourceImage = groundProfile.getImageVertex(0, iSeg);
        Coordinate realReceiver = groundProfile.getVertex(lastVertex);
        Coordinate secondarySource = groundProfile.getVertex(iStart);
        Coordinate secondaryReceiver = groundProfile.getVertex(iEnd);
        Complex[] pImage;
        Complex[] p;
        if (iStart == 0 && iEnd == lastVertex){ // Case 1 no diffraction
            pImage = unitSphericalWavePressure(realSourceImage.distance(realReceiver));
            p = unitSphericalWavePressure(realSource.distance(realReceiver));
        } else {
            p = new Complex[nFreq];
            pImage = new Complex[nFreq];
            if (iStart == 0 && iEnd < lastVertex) { // Case 2 diffraction on the receiver side only
                computeDiffractionAttenuation(realSourceImage, realReceiver, secondaryReceiver, pImage);
                computeDiffractionAttenuation(realSource, realReceiver, secondaryReceiver, p);
            } else if (iStart > 0 && iEnd == lastVertex) { // Case 3 diffraction on the source side only
                Coordinate receiverImage = groundProfile.getImageVertex(lastVertex, iSeg);
                computeDiffractionAttenuation(realSource, receiverImage, secondarySource, pImage);
                computeDiffractionAttenuation(realSource, realReceiver, secondarySource, p);
            } else { // Case 4 diffraction on both sides
                Coordinate secondarySourceImage = groundProfile.getImageVertex(iStart, iSeg);
                Coordinate secondaryReceiverImage = groundProfile.getImageVertex(iEnd, iSeg);
                Complex[] pImage1 = pImage.clone();
                Complex[] p1 = p.clone();
                Complex[] pImage2 = pImage.clone();
                Complex[] p2 = p.clone();
                computeDiffractionAttenuation(realSource, secondaryReceiverImage, secondarySource, pImage1);
                computeDiffractionAttenuation(realSource, secondaryReceiver, secondarySource, p1);
                computeDiffractionAttenuation(secondarySourceImage, realReceiver, secondaryReceiver, pImage2);
                computeDiffractionAttenuation(secondarySource, realReceiver, secondaryReceiver, p2);
                for (int i = 0; i < nFreq; i++) {
                    pImage[i] = pImage1[i].multiply(pImage2[i]);
                    p[i] = p1[i].multiply(p2[i]);
                }
            }
        }
        return (Complex[]) IntStream.range(0, nFreq)
                .mapToObj(i -> pImage[i].divide(p[i]))
                .toArray();
    }

    /**
     * Compute the coherence factor corresponding to a segment
     * Ref: "Coherence factor" subsection of section 2.4.1 from Salomons et al.
     *
     * @param iSeg index of the first index of the segment
     * @param iStart first index of the ground profile or sub-profile
     * @param iEnd last index of the ground profile or sub-profile
     * @return coherence factor
     */
    private double[] coherenceFactor(int iSeg, int iStart, int iEnd){
        Coordinate source = groundProfile.getVertex(iStart);
        Coordinate sourceImage = groundProfile.getImageVertex(iStart, iSeg);
        Coordinate receiver = groundProfile.getVertex(iEnd);
        double sourceHeight = source.getY();
        double receiverHeight = receiver.getY();
        // Standard deviation of the frequency integration
        List<Integer> frequencies = scene.defaultCnossosParameters.getFrequencies();
        double bandwidth = 1./3;
        if (frequencies.get(1) == 2 * frequencies.getFirst()){
            bandwidth = 1;
        }
        double stdFTerm = Math.pow((Math.pow(2, bandwidth/2) - Math.pow(2, -bandwidth/2)) / 3, 2);
        // Standard deviation of the sound speed fluctuations (considered in coherence factor cb below)
        double stdTermC0 = 0;
        // Standard deviation of the source-receiver distance
        double stdDTerm = 0;
        // Standard deviation of the source and receiver heights
        double stdHsTerm = 0;
        if (iStart == 0){ // The source of the portion of profile is the real source
            stdHsTerm = Math.min(1, Math.pow(stdHs / sourceHeight,2));
        }
        double stdHrTerm = 0;
        if (iEnd == groundProfile.getNVertices()-1){ // The source of the portion of profile is the real source
            stdHrTerm = Math.min(1, Math.pow(stdHr / receiverHeight,2));
        }
        // Standard deviation of the phase difference's fluctuation
        double pathLengthDiff = sourceImage.distance(receiver) - source.distance(receiver);
        double stdPhaseTerm = stdFTerm + stdTermC0 + stdDTerm + stdHsTerm + stdHrTerm;
        double[] stdPhase = Arrays.stream(waveNumber).map(k -> Math.sqrt(stdPhaseTerm) * k * pathLengthDiff).toArray();
        // coherence factor Ca
        double[] ca = Arrays.stream(stdPhase).map(s -> Math.exp(-0.5 * s * s )).toArray();
        // coherence factor Cb
        double rho;
        if (sourceHeight == 0 && receiverHeight == 0){
            rho = 0;
        } else {
            rho = sourceHeight * receiverHeight / (sourceHeight + receiverHeight);
        }
        double[] cb = Arrays.stream(waveNumber)
                .map(k -> -3.0/8 * 0.364 * gammaT * k * k * Math.pow(rho, 5.0/3) * source.distance(receiver))
                .toArray();
        return IntStream.range(0, ca.length).mapToDouble(i -> ca[i] * cb[i]).toArray();
    }

    /**
     * Computes the Fresnel ellipse center and semi major axis for a Fresnel parameter of 8
     * Ref: "Fresnel weighting" subsection of section 2.4.2 from Salomons et al.
     *
     * @param iSeg index of the first index of the segment
     * @param iStart first index of the ground profile or sub-profile
     * @param iEnd last index of the ground profile or sub-profile
     * @param center local d coordinate of the ellipse center
     * @param semiMajorAxis ellipse semi major axis length
     */
    private void fresnelEllipse(int iSeg, int iStart, int iEnd, double[] center, double[] semiMajorAxis){
        double[] fresnelParam = new double[waveNumber.length];
        Arrays.fill(fresnelParam, 8);
        fresnelEllipse(iSeg, iStart, iEnd, fresnelParam,center, semiMajorAxis);
    }

    /**
     * Computes the Fresnel ellipse center and semi major axis
     * Ref: "Fresnel weighting" subsection of section 2.4.2 from Salomons et al.
     *
     * @param iSeg index of the first index of the segment
     * @param iStart first index of the ground profile or sub-profile
     * @param iEnd last index of the ground profile or sub-profile
     * @param fresnelParam Fresnel parameter
     * @param center local d coordinate of the ellipse center
     * @param semiMajorAxis ellipse semi major axis length
     */
    private void fresnelEllipse(int iSeg, int iStart, int iEnd, double[] fresnelParam, double[] center, double[] semiMajorAxis){
        Coordinate source = groundProfile.getVertex(iStart);
        Coordinate receiver = groundProfile.getVertex(iEnd);
        Coordinate imageReceiver = groundProfile.getImageVertex(iEnd, iSeg);
        double srcRcvDistance = source.distance(receiver);
        double srcImageReceiverDistance = source.distance(imageReceiver);
        double localSourceHeight = groundProfile.getLocalSourceHeight(iSeg, iStart);
        double localReceiverHeight = groundProfile.getLocalReceiverHeight(iSeg,iEnd );
        double _term = Math.sqrt(Math.pow(localSourceHeight + localReceiverHeight, 2) + Math.pow(srcImageReceiverDistance, 2));
        double[] d = IntStream.range(0, waveNumber.length)
                .mapToDouble(i -> 2 * Math.PI / waveNumber[i] / fresnelParam[i] + _term)
                .toArray(); // Eq. 41
        center = Arrays.stream(d)
                .map(di -> srcRcvDistance / 2 * (Math.pow(localSourceHeight,2) - Math.pow(localReceiverHeight,2))
                        / (Math.pow(di,2) - Math.pow(srcRcvDistance,2)))
                .toArray();
        double[] dsSquare = Arrays.stream(center)
                .map(df -> Math.pow(df, 2) + Math.pow(localSourceHeight,2))
                .toArray();
        double[] drSquare = Arrays.stream(center)
                .map(df -> Math.pow(srcRcvDistance - df, 2) + Math.pow(localReceiverHeight,2))
                .toArray();
        semiMajorAxis = IntStream.range(0, d.length)
                .mapToDouble(i -> 0.5 * Math.sqrt(
                        (Math.pow(d[i],4) + Math.pow(dsSquare[i] - drSquare[i], 2) - 2 * Math.pow(d[i],2) * (dsSquare[i] + drSquare[i]))
                        / (Math.pow(d[i],2) - Math.pow(srcRcvDistance,2))
                ))
                .toArray();
    }
}
