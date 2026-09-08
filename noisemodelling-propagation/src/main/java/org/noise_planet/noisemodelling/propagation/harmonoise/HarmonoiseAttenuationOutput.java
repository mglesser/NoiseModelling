/**
 * NoiseModelling is a library capable of producing noise maps. It can be freely used either for research and education, as well as by experts in a professional use.
 * <p>
 * NoiseModelling is distributed under GPL 3 license. You can read a copy of this License in the file LICENCE provided with this software.
 * <p>
 * Official webpage : http://noise-planet.org/noisemodelling.html
 * Contact: contact@noise-planet.org
 */

package org.noise_planet.noisemodelling.propagation.harmonoise;

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
    HarmonoiseGroundProfile groundProfile;
    double excessAttenuation = 0;
    List<Double> diffractionAttenuation = new ArrayList<>();

    public HarmonoiseAttenuationOutput(HarmonoiseGroundProfile groundProfile){
        this.groundProfile = groundProfile;
    }

    /**
     * Add frequency dependant attenuation to diffractionAttenuation
     *
     * @param attenuation attenuation to add
     */
    public void addDiffractionAttenuation(List<Double> attenuation){
        if (diffractionAttenuation.isEmpty()) {
            diffractionAttenuation = attenuation;
        } else {
            diffractionAttenuation = IntStream.range(0, diffractionAttenuation.size())
                    .mapToObj(i -> diffractionAttenuation.get(i) + attenuation.get(i))
                    .collect(Collectors.toList());
        }
    }
}
