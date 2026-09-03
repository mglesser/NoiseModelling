/**
 * NoiseModelling is a library capable of producing noise maps. It can be freely used either for research and education, as well as by experts in a professional use.
 * <p>
 * NoiseModelling is distributed under GPL 3 license. You can read a copy of this License in the file LICENCE provided with this software.
 * <p>
 * Official webpage : http://noise-planet.org/noisemodelling.html
 * Contact: contact@noise-planet.org
 */

package org.noise_planet.noisemodelling.propagation.harmonoise;

import org.noise_planet.noisemodelling.propagation.AttenuationParameters;
import org.noise_planet.noisemodelling.propagation.SceneWithAttenuation;
import org.noise_planet.noisemodelling.propagation.cnossos.CnossosAttenuationOutput;

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
        computeExcessAttenuation(data, scene, attenuationOutput, exportAttenuationMatrix, 0,
                attenuationOutput.groundProfile.profile.getNumPoints());
    }

    /**
     * Recursive calculation scheme for excess attenuation
     * Ref: section 2.2.3 from Salomons et al.
     *
     * @param data Attenuation parameters
     * @param scene Scene with attenuation data
     * @param attenuationOutput Output of the attenuation computation
     * @param exportAttenuationMatrix if true, store intermediate values in attenuationOutput for debugging purpose
     * @param startIndex index of the first ground profile vertex to consider
     * @param endIndex index of the last ground profile vertex to consider
     */
    public static void computeExcessAttenuation(AttenuationParameters data, SceneWithAttenuation scene,
                                                HarmonoiseAttenuationOutput attenuationOutput,
                                                boolean exportAttenuationMatrix, int startIndex, int endIndex) {

    }
}
