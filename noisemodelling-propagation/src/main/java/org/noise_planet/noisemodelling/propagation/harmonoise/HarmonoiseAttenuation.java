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

/**
 * Compute excess attenuation according to Harmonoise propagation model
 * Ref: Salomons, E., Van Maercke, D., Defrance, J.,&amp;De Roo, F. (2011). The Harmonoise sound propagation model.
 * Acta acustica united with acustica, 97(1), 62-74
 * @author Martin Glesser
 */

public class HarmonoiseAttenuation {
    SceneWithAttenuation scene; // Scene with attenuation data
    HarmonoiseAttenuationOutput attenuationOutput; // Output of the attenuation computation
    //    public boolean exportAttenuationMatrix; // if true, store intermediate values for debugging purpose

    public HarmonoiseAttenuation(SceneWithAttenuation scene, HarmonoiseAttenuationOutput output) {
        this.scene = scene;
        this.attenuationOutput = output;
    }

    /**
     * Recursive calculation scheme for excess attenuation
     * Ref: section 2.2.3 from Salomons et al.
     */
    public void computeExcessAttenuation() {
        Coordinate[] vertices = attenuationOutput.groundProfile.vertices;
        vertices[0].y = attenuationOutput.groundProfile.getSource().z;
        vertices[vertices.length-1].y = attenuationOutput.groundProfile.getReceiver().z;
        computeExcessAttenuation(vertices);
    }

    /**
     * Recursive calculation scheme for excess attenuation
     * Ref: section 2.2.3 from Salomons et al.
     *
     * @param vertices ground vertices (incl. src and rcv)
     */
    private void computeExcessAttenuation(Coordinate[] vertices){
        double maxDistance = 0;
        int indexMaxDistance = 0;
        for (int i = 1; i < vertices.length - 1; i++) {
            double crossProduct = (vertices[i].x - vertices[0].x) * (vertices[0].y - vertices[vertices.length-1].y)
                    + (vertices[i].y - vertices[0].y) * (vertices[vertices.length-1].x - vertices[0].x);
            if (crossProduct > 0) { // if the current point is above the line crossing the first and last points
                double distance = vertices[0].distance(vertices[2]) + vertices[2].distance(vertices[vertices.length-1])
                        - vertices[0].distance(vertices[vertices.length-1]); // Eq. 5
                if (distance > maxDistance){
                    maxDistance = distance;
                    indexMaxDistance = i;
                }
            }
        }
        if (indexMaxDistance == 0) { // No diffraction point above the line crossing the first and last points
            computeGroundAttenuation(vertices);
        }else {
            computeDiffractionAttenuation(vertices[0], vertices[vertices.length-1],
                    vertices[indexMaxDistance]);
            computeExcessAttenuation(Arrays.copyOfRange(vertices, 0, indexMaxDistance+1));
            computeExcessAttenuation(Arrays.copyOfRange(vertices, indexMaxDistance, vertices.length));
        }
    }

    /**
     * Compute diffraction attenuation according to Harmonoise methodology
     * Ref: section 2.3 from Salomons et al.
     *
     * @param source source point ("real" or secondary at diffraction edge)
     * @param receiver receiver point ("real" or secondary at diffraction edge)
     * @param point diffraction point
     */
    private void computeDiffractionAttenuation(Coordinate source, Coordinate receiver, Coordinate point){
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
        List<Double> fresnelNumber = scene.defaultCnossosParameters.getFrequenciesExact()
                .stream()
                .map(f -> 2 * pathLengthDiff / (scene.defaultCnossosParameters.getCelerity() / f))
                .toList(); // Eq.8
        List<Double> diffractionAttenuation = fresnelNumber.stream()
                .map(HarmonoiseAttenuation::fresnelApproximation)
                .toList();
        attenuationOutput.addDiffractionAttenuation(diffractionAttenuation);
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

    private void computeGroundAttenuation(Coordinate[] vertices){
        if (hasConvexSegment(vertices)){
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
}
