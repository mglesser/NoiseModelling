/**
 * NoiseModelling is a library capable of producing noise maps. It can be freely used either for research and education, as well as by experts in a professional use.
 * <p>
 * NoiseModelling is distributed under GPL 3 license. You can read a copy of this License in the file LICENCE provided with this software.
 * <p>
 * Official webpage : http://noise-planet.org/noisemodelling.html
 * Contact: contact@noise-planet.org
 */

package org.noise_planet.noisemodelling.propagation.harmonoise;

import org.noise_planet.noisemodelling.pathfinder.profilebuilder.CutProfile;
import org.noise_planet.noisemodelling.propagation.AttenuationOutput;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Output of the Harmonoise attenuation computation (data class).
 *
 * @author Martin Glesser
 */
public class HarmonoiseAttenuationOutput extends AttenuationOutput {
    double excessAttenuation = 0;
    double[] diffractionAttenuation = new double[0];

    public HarmonoiseAttenuationOutput(CutProfile cutProfile){
        this.cutProfile = cutProfile;
    }

    /**
     * Add frequency dependant attenuation to diffractionAttenuation
     *
     * @param attenuation attenuation to add
     */
    public void addDiffractionAttenuation(double[] attenuation){
        if (diffractionAttenuation.length == 0) {
            diffractionAttenuation = attenuation;
        } else {
            diffractionAttenuation = IntStream.range(0, diffractionAttenuation.length)
                    .mapToDouble(i -> diffractionAttenuation[i] + attenuation[i])
                    .toArray();
        }
    }
}
