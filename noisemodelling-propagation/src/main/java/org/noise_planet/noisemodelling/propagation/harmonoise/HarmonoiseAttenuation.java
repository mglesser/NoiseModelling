/**
 * NoiseModelling is a library capable of producing noise maps. It can be freely used either for research and education, as well as by experts in a professional use.
 * <p>
 * NoiseModelling is distributed under GPL 3 license. You can read a copy of this License in the file LICENCE provided with this software.
 * <p>
 * Official webpage : http://noise-planet.org/noisemodelling.html
 * Contact: contact@noise-planet.org
 */

package org.noise_planet.noisemodelling.propagation.harmonoise;

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

    /**
     * Recursive calculation scheme for excess attenuation
     * Ref: section 2.2.3 from Salomons et al.
     *
     * @param scene Scene with attenuation data
     * @param attenuationOutput Output of the attenuation computation
     * @param exportAttenuationMatrix if true, store intermediate values in attenuationOutput for debugging purpose
     */
    public static void computeExcessAttenuation(SceneWithAttenuation scene,
                                                 HarmonoiseAttenuationOutput attenuationOutput,
                                                boolean exportAttenuationMatrix) {
        Coordinate[] vertices = attenuationOutput.groundProfile.vertices;
        vertices[0].y = attenuationOutput.groundProfile.getSource().z;
        vertices[vertices.length-1].y = attenuationOutput.groundProfile.getReceiver().z;
        computeExcessAttenuation(scene, attenuationOutput, vertices);
    }

    /**
     * Recursive calculation scheme for excess attenuation
     * Ref: section 2.2.3 from Salomons et al.
     *
     * @param scene Scene with attenuation data
     * @param vertices ground vertices (incl. src and rcv)
     */
    public static void computeExcessAttenuation(SceneWithAttenuation scene,
                                                HarmonoiseAttenuationOutput attenuationOutput, Coordinate[] vertices){
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
            computeGroundAttenuation(attenuationOutput);
        }else {
            computeDiffractionAttenuation(scene, attenuationOutput, vertices[0], vertices[vertices.length-1],
                    vertices[indexMaxDistance]);
            computeExcessAttenuation(scene, attenuationOutput,
                    Arrays.copyOfRange(vertices, 0, indexMaxDistance+1));
            computeExcessAttenuation(scene, attenuationOutput,
                    Arrays.copyOfRange(vertices, indexMaxDistance, vertices.length));
        }
    }

    public static void computeGroundAttenuation(HarmonoiseAttenuationOutput attenuationOutput){
        attenuationOutput.excessAttenuation += 0;
    }

    /**
     * Compute diffraction attenuation according to Harmonoise methodology
     * Ref: section 2.3 from Salomons et al.
     *
     * @param scene Scene with attenuation data
     * @param attenuationOutput Output of the attenuation computation
     * @param source source point ("real" or intermediate)
     * @param receiver receiver point ("real" or intermediate)
     * @param point diffraction point
     */
    public static void computeDiffractionAttenuation(SceneWithAttenuation scene,
                                                     HarmonoiseAttenuationOutput attenuationOutput, Coordinate source,
                                                     Coordinate receiver, Coordinate point){
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
    public static double fresnelApproximation(double fresnelNumber){
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

}
