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
    private final double[] waveNumber;

    // To be exposed to user (input of the model)
    double stdHs = 0; // standard deviation on the source height
    double stdHr = 0; // standard deviation on the receiver height
    double gammaT = 5 * Math.pow(10,-6); // typical value for moderate turbulence

    public HarmonoiseAttenuation(SceneWithAttenuation scene, HarmonoiseAttenuationOutput output) {
        this.scene = scene;
        this.attenuationOutput = output;
        waveNumber = scene.defaultCnossosParameters.getFrequenciesExact()
                .stream()
                .mapToDouble(f -> 2 * Math.PI * f / scene.defaultCnossosParameters.getCelerity())
                .toArray();
        this.groundProfile = new HarmonoiseGroundProfile(output.getCutProfile(), waveNumber.length);
    }

    /**
     * Recursive calculation scheme for excess attenuation
     * Ref: section 2.2.3 from Salomons et al.
     */
    public void computeExcessAttenuation() {
        computeExcessAttenuation(0, groundProfile.getNVertices()-1);
    }

    /**
     * Compute the excess attenuation for a ground profile from (secondary) source to (secondary) receiver
     * Ref: section 2.2.3 from Salomons et al.
     *
     * @param iSource index of the (secondary) source
     * @param iReceiver index of the (secondary) receiver
     */
    private void computeExcessAttenuation(int iSource, int iReceiver){
        Coordinate[] vertices = groundProfile.getVertices(iSource, iReceiver);
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
            computeGroundAttenuation(iSource, iReceiver);
        }else {
            computeDiffractionAttenuation(vertices[0], vertices[vertices.length-1],
                    vertices[indexMaxDistance]);
            computeExcessAttenuation(iSource, iSource + indexMaxDistance);
            computeExcessAttenuation(iSource + indexMaxDistance, iReceiver);
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

    private void computeGroundAttenuation(int iSource, int iReceiver) {

        if (hasConvexSegment(groundProfile.getVertices(iSource, iReceiver))){
            attenuationOutput.excessAttenuation += 0;
        } else {
            concaveGroundAttenuation(iSource, iReceiver);
        }

    }

    private void concaveGroundAttenuation(int iSource, int iReceiver){
        double transitionFrequency = transitionFrequency(iSource, iReceiver);
        // Pre-computation
        Complex[][] geometricalWeightingFactor = new Complex[groundProfile.getNVertices()-1][waveNumber.length];
        double[][] coherenceFactor = new double[groundProfile.getNVertices()-1][waveNumber.length];
        double[][] modifiedFresnelWeighting = new double[groundProfile.getNVertices()-1][waveNumber.length];
        for (int k = iSource; k < iReceiver; k++) {
            groundProfile.setReflectionCoefficient(k, sphericalWaveReflectionCoefficient(k, iSource, iReceiver));
            geometricalWeightingFactor[k] = geometricalWeightingFactor(k, iSource, iReceiver);
            coherenceFactor[k] = coherenceFactor(k, iSource, iReceiver);
            modifiedFresnelWeighting[k] = modifiedFresnelWeighting(k, iSource, iReceiver, transitionFrequency);
        }

    }

    /**
     * Return the ground attenuation of a (sub)profile for relatively flat ground
     * Ref: Eq. 19 and 20 from Salomons et al.
     *
     * @param geometricalWeightingFactor geometrical weighting factors for all segments of the (sub)profile
     * @param coherenceFactor coherence factor for all segments of the (sub)profile
     * @param modifiedFresnelWeighting modified Fresnel weighing for all segments of the (sub)profile
     * @return round attenuation for relatively flat ground
     */
    private double[] flatGroundAttenuation(Complex[][] geometricalWeightingFactor, double[][] coherenceFactor,
                                           double[][] modifiedFresnelWeighting) {
        int nbSeg = geometricalWeightingFactor.length;
        int nbFreq = geometricalWeightingFactor[0].length;
        double[] flatGroundAttenuation = new double[nbFreq];
        Arrays.fill(flatGroundAttenuation, 0);
        for (int iSeg = 0; iSeg < nbSeg; iSeg++) {
            int finalISeg = iSeg;
            double[] attenuation = IntStream.range(0, nbFreq)
                    .mapToDouble(i -> flatGroundAttenuation[i]
                            + modifiedFresnelWeighting[finalISeg][i] * 10 * Math.log10(
                                    Math.pow(
                                            groundProfile.getReflectionCoefficient(finalISeg, i)
                                                    .multiply(geometricalWeightingFactor[finalISeg][i])
                                                    .multiply(coherenceFactor[finalISeg][i])
                                                    .add(1).abs(),2
                                    ) + Math.pow(
                                            groundProfile.getReflectionCoefficient(finalISeg, i)
                                                    .multiply(geometricalWeightingFactor[finalISeg][i])
                                                    .abs(),2
                                    ) * (1 - Math.pow(coherenceFactor[finalISeg][i],2))
                            )
                    )
                    .toArray();
        }
        return flatGroundAttenuation;
    }

    /**
     * Determine whether the ground profile contains convex segment or not
     * Ref: Figure 3 and "Geometry" subsection of section 2.4.1 from Salomons et al.
     *
     * @param vertices vertices of the ground profile
     * @return true if the profile contains at least one convex segment
     */
    private static boolean hasConvexSegment(Coordinate[] vertices){
        int iReceiver = vertices.length - 1;
        boolean isConvex = false;
        // Loop on segments
        for (int i = 0; i < vertices.length - 2; i++) {
            double localSourceHeight = vertices[i+1].distance(vertices[0])
                    * Math.sin(Angle.angleBetweenOriented(vertices[0], vertices[i+1], vertices[i]));
            double localReceiverHeight = vertices[i].distance(vertices[iReceiver])
                    * Math.sin(Angle.angleBetweenOriented(vertices[i+1], vertices[i], vertices[iReceiver]));
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
     * @param iSeg index of the first index of the reflection plane segment
     * @param iSource index of the (secondary) source
     * @param iReceiver index of the (secondary) receiver
     * @return spherical wave reflection coefficient
     */
    private Complex[] sphericalWaveReflectionCoefficient(int iSeg, int iSource, int iReceiver){
        double angle = groundProfile.getReflexionAngle(iSeg, iSource, iReceiver);
        List<Double> frequencies = scene.defaultCnossosParameters.getFrequenciesExact();
        double distance = groundProfile.getImageVertex(iSource, iSeg).distance(groundProfile.getVertex(iReceiver));
        double hm = (groundProfile.getLocalOrdinate(iSource, iSeg)
                + groundProfile.getLocalOrdinate(iReceiver, iSeg)) / 2;
        Complex[] groundImpedance = groundProfile.getGroundImpedance(iSeg, frequencies);
        Complex[] reflectionCoefficient = new Complex[frequencies.size()];
        for (int i = 0; i < frequencies.size(); i++) {
            Complex z = groundImpedance[i].multiply(Math.cos(angle));
            Complex planeWaveReflectionCoefficient = z.subtract(1).divide(z.add(1));
            double kr = 2 * Math.PI * frequencies.get(i) / scene.defaultCnossosParameters.getCelerity() * distance;
            Complex admittance = groundImpedance[i].reciprocal();
            Complex numericalDistance = admittance.multiply(Math.sqrt(kr))
                    .multiply(new Complex(0.5, 0.5));
            double hG = scene.defaultCnossosParameters.celerity / frequencies.get(i) / 32;
            double nG = 1 - 0.7 * Math.exp(-hm / hG);
            reflectionCoefficient[i] = boundaryLossFactor(numericalDistance).pow(nG)
                    .multiply(new Complex(1).subtract(planeWaveReflectionCoefficient))
                    .add(planeWaveReflectionCoefficient);
        }
        return reflectionCoefficient;
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
     * @param iSource index of the (secondary) source
     * @param iReceiver index of the (secondary) receiver
     * @return geometrical weighting factor
     */
    private Complex[] geometricalWeightingFactor(int iSeg, int iSource, int iReceiver){
        int nFreq = scene.defaultCnossosParameters.getFrequenciesExact().size();
        int lastVertex = groundProfile.getNVertices()-1;
        Coordinate realSource = groundProfile.getVertex(0);
        Coordinate realSourceImage = groundProfile.getImageVertex(0, iSeg);
        Coordinate realReceiver = groundProfile.getVertex(lastVertex);
        Coordinate secondarySource = groundProfile.getVertex(iSource);
        Coordinate secondaryReceiver = groundProfile.getVertex(iReceiver);
        Complex[] pImage;
        Complex[] p;
        if (iSource == 0 && iReceiver == lastVertex){ // Case 1 no diffraction
            pImage = unitSphericalWavePressure(realSourceImage.distance(realReceiver));
            p = unitSphericalWavePressure(realSource.distance(realReceiver));
        } else {
            p = new Complex[nFreq];
            pImage = new Complex[nFreq];
            if (iSource == 0 && iReceiver < lastVertex) { // Case 2 diffraction on the receiver side only
                computeDiffractionAttenuation(realSourceImage, realReceiver, secondaryReceiver, pImage);
                computeDiffractionAttenuation(realSource, realReceiver, secondaryReceiver, p);
            } else if (iSource > 0 && iReceiver == lastVertex) { // Case 3 diffraction on the source side only
                Coordinate receiverImage = groundProfile.getImageVertex(lastVertex, iSeg);
                computeDiffractionAttenuation(realSource, receiverImage, secondarySource, pImage);
                computeDiffractionAttenuation(realSource, realReceiver, secondarySource, p);
            } else { // Case 4 diffraction on both sides
                Coordinate secondarySourceImage = groundProfile.getImageVertex(iSource, iSeg);
                Coordinate secondaryReceiverImage = groundProfile.getImageVertex(iReceiver, iSeg);
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
     * @param iSource index of the (secondary) source
     * @param iReceiver index of the (secondary) receiver
     * @return coherence factor
     */
    private double[] coherenceFactor(int iSeg, int iSource, int iReceiver){
        Coordinate source = groundProfile.getVertex(iSource);
        Coordinate sourceImage = groundProfile.getImageVertex(iSource, iSeg);
        Coordinate receiver = groundProfile.getVertex(iReceiver);
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
        if (iSource == 0){ // The source of the portion of profile is the real source
            stdHsTerm = Math.min(1, Math.pow(stdHs / sourceHeight,2));
        }
        double stdHrTerm = 0;
        if (iReceiver == groundProfile.getNVertices()-1){ // The source of the portion of profile is the real source
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
     * @param iSource index of the (secondary) source
     * @param iReceiver index of the (secondary) receiver
     * @param iFreq index of the frequency
     * @param center local d coordinate of the ellipse center
     * @param semiMajorAxis ellipse semi major axis length
     */
    private void fresnelEllipse(int iSeg, int iSource, int iReceiver, int iFreq, double center, double semiMajorAxis){
        fresnelEllipse(iSeg, iSource, iReceiver, iFreq, 8,center, semiMajorAxis);
    }

    /**
     * Computes the Fresnel ellipse center and semi major axis
     * Ref: "Fresnel weighting" subsection of section 2.4.2 from Salomons et al.
     *
     * @param iSeg index of the first index of the segment
     * @param iSource index of the (secondary) source
     * @param iReceiver index of the (secondary) receiver
     * @param iFreq index of the frequency
     * @param fresnelParam Fresnel parameter
     * @param center local d coordinate of the ellipse center
     * @param semiMajorAxis ellipse semi major axis length
     */
    private void fresnelEllipse(int iSeg, int iSource, int iReceiver, int iFreq, double fresnelParam, double center, double semiMajorAxis){
        Coordinate source = groundProfile.getVertex(iSource);
        Coordinate receiver = groundProfile.getVertex(iReceiver);
        Coordinate imageReceiver = groundProfile.getImageVertex(iReceiver, iSeg);
        double srcRcvDistance = source.distance(receiver);
        double srcImageReceiverDistance = source.distance(imageReceiver);
        double localSourceHeight = groundProfile.getLocalOrdinate(iSource, iSeg);
        double localReceiverHeight = groundProfile.getLocalOrdinate(iReceiver, iSeg);
        double _term = Math.sqrt(Math.pow(localSourceHeight + localReceiverHeight, 2) + Math.pow(srcImageReceiverDistance, 2));
        double d = 2 * Math.PI / waveNumber[iFreq] / fresnelParam + _term; // Eq. 41
        double denominator = Math.pow(d, 2) - Math.pow(srcRcvDistance, 2);
        center = srcRcvDistance / 2 * (Math.pow(localSourceHeight,2) - Math.pow(localReceiverHeight,2))
                / denominator;
        double dsSquare = Math.pow(center, 2) + Math.pow(localSourceHeight,2);
        double drSquare = Math.pow(srcRcvDistance - center, 2) + Math.pow(localReceiverHeight,2);
        semiMajorAxis = 0.5 * Math.sqrt(
                (Math.pow(d,4) + Math.pow(dsSquare - drSquare, 2) - 2 * Math.pow(d,2) * (dsSquare + drSquare))
                        / denominator);
    }

    /**
     * Fresnel weighting function used for Fresnel weighting
     * Ref: Eq. 45 from Salomons et al.
     *
     * @param x input parameter
     * @return Fresnel weighting
     */
    private double fresnelFunction(double x){
        double output;
            if (x <= 1){
                output = 0;
            } else if (x >= 1){
                output = 1;
            } else {
                output = 1 - 1/Math.PI * (Math.acos(x) - x * Math.sqrt(1 - Math.pow(x, 2)));
            }
        return output;
    }

    /**
     * Computes the Fresnel weightings of a ground segment at a given frequency
     * Ref: "Fresnel weighting" subsection of section 2.4.2 from Salomons et al.
     *
     * @param iSeg index of the first index of the segment
     * @param iSource index of the (secondary) source
     * @param iReceiver index of the (secondary) receiver
     * @param xi1 input parameter of the Fresnel function
     * @param xi2 input parameter of the Fresnel function
     */
    private double fresnelWeighting(int iSeg, int iSource, int iReceiver, double xi1, double xi2){
        double f1;
        double f2;
        if (iSeg == iSource){
            f1 = 0;
        } else {
            f1 = fresnelFunction(xi1);
        }
        if (iSeg+1 == iReceiver){
            f2 = 1;
        } else {
            f2 = fresnelFunction(xi2);
        }
        return f2 - f1;
    }

    /**
     * Return the frequency dependant modified Fresnel weightings for a given (sub)profile
     * Ref: "Modified Fresnel weighting" subsection of section 2.4.2 from Salomons et al.
     *
     * @param iSeg index of the ground segment
     * @param iSource index of the (secondary) source
     * @param iReceiver index of the (secondary) receiver
     * @param transitionFrequency transition frequency
     * @return modified Fresnel weightings
     */
    private double[] modifiedFresnelWeighting(int iSeg, int iSource, int iReceiver, double transitionFrequency){
        Coordinate source = groundProfile.getVertex(iSource);
        Coordinate receiver = groundProfile.getVertex(iReceiver);
        double localSourceHeight = groundProfile.getLocalOrdinate(iSource, iSeg);
        double localReceiverHeight = groundProfile.getLocalOrdinate(iReceiver, iSeg);
        double dsr = source.distance(receiver);
        double dsp = dsr * localReceiverHeight / (localSourceHeight + localReceiverHeight);
        List<Double> frequency = scene.defaultCnossosParameters.getFrequenciesExact();
        double[] nf = frequency.stream()
                .mapToDouble(f -> 32 * (1 - Math.exp(Math.pow(transitionFrequency,2)/Math.pow(f,2))))
                .toArray();
        double[] modifiedWeighting = new double[frequency.size()];
        double center = 1;
        double semiMajorAxis = 1;
        double alpha, dc, xiC, xi1, xi2, xi1Prim, xi2Prim;
        for (int i = 0; i < frequency.size(); i++) {
            fresnelEllipse(iSeg, iSource, iReceiver, i, nf[i], center, semiMajorAxis);
            alpha = Math.pow((1 + Math.pow(frequency.get(i) / transitionFrequency, 2)), -1);
            dc = alpha * center + (1 - alpha) * dsp;
            xiC = (dc - center) / semiMajorAxis;
            xi1 = (groundProfile.getLocalAbscissa(iSeg, iSeg, iSource) - center) / semiMajorAxis;
            xi2 = (groundProfile.getLocalAbscissa(iSeg+1, iSeg, iSource) - center) / semiMajorAxis;
            xi1Prim = (xi1 - xiC) / (1 - xi1*xiC);
            xi2Prim = (xi2 - xiC) / (1 - xi2*xiC);
            modifiedWeighting[i] = fresnelWeighting(iSeg, iSource,iReceiver, xi1Prim, xi2Prim);
        }
    return modifiedWeighting;
    }

    /**
     *  Return, for a given (sub)profile and reflection plane, the phase difference between direct and reflected sound
     * /!\ groundProfile.reflectionCoefficient must have been pre-computed before calling the present method
     * Ref: Eq. 57 from Salomons et al.
     *
     * @param iSeg index of the reflection plane segment
     * @param iSource index of the (secondary) source
     * @param iReceiver index of the (secondary) receiver
     * @param iFreq index of the frequency
     * @return phase difference
     */
    private double phaseDifference(int iSeg, int iSource, int iReceiver, int iFreq){
        Coordinate source = groundProfile.getVertex(iSource);
        Coordinate receiver = groundProfile.getVertex(iReceiver);
        Coordinate imageSource = groundProfile.getImageVertex(iSource, iSeg);
        Complex reflectionCoefficient = groundProfile.getReflectionCoefficient(iSeg, iFreq);
        return reflectionCoefficient.getArgument() +
                waveNumber[iFreq] * (imageSource.distance(receiver) - source.distance(receiver));
    }

    /**
     * Return the transition frequency for a given (sub)profile
     * /!\ groundProfile.reflectionCoefficient must have been pre-computed before calling the present method
     * Eq. 53 to 58 from Salomons et al.
     *
     * @param iSource index of the (secondary) source
     * @param iReceiver index of the (secondary) receiver
     * @return transition frequency
     */
    double transitionFrequency(int iSource, int iReceiver){
        List<Double> frequency = scene.defaultCnossosParameters.getFrequenciesExact();
        double fMin = 0;
        double fMax = 0;
        double maxPhaseDiff = 0;
        double phaseDiffPrevious = 0;
        double phaseDiff = 0;
        double center = 1;
        double semiMajorAxis = 1;
        for (int i = 0; i < frequency.size(); i++) {
            for (int k = iSource; k < iReceiver; k++) {
                if (groundProfile.isConvexSegment(k, iSource, iReceiver)){
                    continue;
                }
                fresnelEllipse(k, iSource, iReceiver, i, center, semiMajorAxis);
                double xi1 = (groundProfile.getLocalAbscissa(k, k, iSource) - center) / semiMajorAxis;
                double xi2 = (groundProfile.getLocalAbscissa(k+1, k, iSource) - center) / semiMajorAxis;
                if (fresnelWeighting(k, iSource, iReceiver,xi1, xi2) == 0){
                    continue;
                }
                phaseDiffPrevious = phaseDiff;
                phaseDiff = phaseDifference(k, iSource, iReceiver, i);
                if (phaseDiff > maxPhaseDiff){
                    maxPhaseDiff = phaseDiff;
                }
            }
            if (maxPhaseDiff > Math.PI/2 && fMin == 0) {
                if (i == 0){
                    fMin = frequency.getFirst() / 2;
                } else if (i == frequency.size()){
                    fMin = frequency.getLast() ;
                } else {
                    fMin = frequency.get(i-1) + (frequency.get(i) - frequency.get(i-1))
                            * (Math.PI/2 - phaseDiffPrevious) / (phaseDiff - phaseDiffPrevious) ;
                }
            }
            if (maxPhaseDiff > Math.PI) {
                if (i == 0){
                    fMin = frequency.getFirst();
                } else if (i == frequency.size()){
                    fMin = frequency.get(frequency.size()-2) * 2 ;
                } else {
                    fMin = frequency.get(i-1) + (frequency.get(i) - frequency.get(i-1))
                            * (Math.PI - phaseDiffPrevious) / (phaseDiff - phaseDiffPrevious) ;
                }
                break;
            }
        }
        return Math.sqrt(fMin * fMax);
    }
}
