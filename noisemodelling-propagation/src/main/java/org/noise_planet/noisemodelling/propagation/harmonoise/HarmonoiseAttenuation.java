/**
 * NoiseModelling is a library capable of producing noise maps. It can be freely used either for research and education, as well as by experts in a professional use.
 * <p>
 * NoiseModelling is distributed under GPL 3 license. You can read a copy of this License in the file LICENCE provided with this software.
 * <p>
 * Official webpage : http://noise-planet.org/noisemodelling.html
 * Contact: contact@noise-planet.org
 */

package org.noise_planet.noisemodelling.propagation.harmonoise;

import org.locationtech.jts.geom.Coordinate;
import org.noise_planet.noisemodelling.propagation.AttenuationParameters;
import org.noise_planet.noisemodelling.propagation.SceneWithAttenuation;

import java.util.Arrays;

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
     * @param data Attenuation parameters
     * @param scene Scene with attenuation data
     * @param attenuationOutput Output of the attenuation computation
     * @param exportAttenuationMatrix if true, store intermediate values in attenuationOutput for debugging purpose
     */
    public static void computeExcessAttenuation(AttenuationParameters data, SceneWithAttenuation scene,
                                                 HarmonoiseAttenuationOutput attenuationOutput,
                                                boolean exportAttenuationMatrix) {
        Coordinate[] vertices = attenuationOutput.groundProfile.vertices;
        vertices[0].y = attenuationOutput.groundProfile.getSource().z;
        vertices[vertices.length-1].y = attenuationOutput.groundProfile.getReceiver().z;
        computeExcessAttenuation(data, scene, attenuationOutput, vertices);
    }

    /**
     * Recursive calculation scheme for excess attenuation
     * Ref: section 2.2.3 from Salomons et al.
     *
     * @param data Attenuation parameters
     * @param scene Scene with attenuation data
     * @param vertices ground vertices (incl. src and rcv)
     */
    public static void computeExcessAttenuation(AttenuationParameters data, SceneWithAttenuation scene,
                                                HarmonoiseAttenuationOutput attenuationOutput, Coordinate[] vertices){
        double maxDistance = 0;
        int indexMaxDistance = 0;
        for (int i = 1; i < vertices.length - 1; i++) {
            double crossProduct = (vertices[i].x - vertices[0].x) * (vertices[0].y - vertices[vertices.length-1].y)
                    + (vertices[i].y - vertices[0].y) * (vertices[vertices.length-1].x - vertices[0].x);
            if (crossProduct > 0) { // if the current point is above the line crossing the first and last points
                double distance = vertices[0].distance(vertices[2]) + vertices[2].distance(vertices[vertices.length-1])
                        - vertices[0].distance(vertices[vertices.length-1]);
                if (distance > maxDistance){
                    maxDistance = distance;
                    indexMaxDistance = i;
                }
            }
        }
        if (indexMaxDistance == 0) { // No diffraction point above the line crossing the first and last points
            computeGroundAttenuation(attenuationOutput);
        }else {
            computeDiffractionAttenuation(attenuationOutput);
            computeExcessAttenuation(data, scene, attenuationOutput,
                    Arrays.copyOfRange(vertices, 0, indexMaxDistance+1));
            computeExcessAttenuation(data, scene, attenuationOutput,
                    Arrays.copyOfRange(vertices, indexMaxDistance, vertices.length));
        }
    }

    public static void computeGroundAttenuation(HarmonoiseAttenuationOutput attenuationOutput){
        attenuationOutput.excessAttenuation += 0;
    }

    public static void computeDiffractionAttenuation(HarmonoiseAttenuationOutput attenuationOutput){
        attenuationOutput.excessAttenuation += 0;
    }

}
