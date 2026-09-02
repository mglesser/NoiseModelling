package org.noise_planet.noisemodelling.propagation;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.junit.jupiter.api.Test;
import org.noise_planet.noisemodelling.pathfinder.profilebuilder.CutProfile;
import org.noise_planet.noisemodelling.pathfinder.profilebuilder.ProfileBuilder;
import org.noise_planet.noisemodelling.propagation.harmonoise.HarmonoisePropagationModel;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class AttenuationComputeOutputHarmonoiseTest {
    private static final double HUMIDITY = 70;
    private static final double TEMPERATURE = 10;

    private static CutProfile loadCutProfile(InputStream inputStream) throws IOException {
        ObjectMapper mapper = new ObjectMapper();
        return mapper.readValue(inputStream, CutProfile.class);
    }

    private static double[] computeHarmonoiseAttenuation(String utName)
            throws IOException {
        //Get test data
        URL url = AttenuationComputeOutputHarmonoiseTest.class.getResource("harmonoise/" + utName + ".json");

        //Create profile builder
        ProfileBuilder profileBuilder = new ProfileBuilder()
                .finishFeeding();

        //Propagation data building
        SceneWithAttenuation sceneWithAttenuation = new SceneWithAttenuation(profileBuilder);
        sceneWithAttenuation.sourceGs.put(-1L, 0.5);

        //Propagation process path data building
        sceneWithAttenuation.defaultCnossosParameters.setHumidity(HUMIDITY);
        sceneWithAttenuation.defaultCnossosParameters.setTemperature(TEMPERATURE);

        //Out and computation settings
        CutProfile cutProfile;
        try(InputStream inputStream = url.openStream()) {
            cutProfile = loadCutProfile(inputStream);
        }

        PropagationModel propagationModel = new HarmonoisePropagationModel();

        return propagationModel.computeAttenuation(sceneWithAttenuation, cutProfile,
                sceneWithAttenuation.defaultCnossosParameters, false).getFirst().getaGlobal();
    }

    /**
     * Test case 1 from Harmonoise publication (hard/rigid ground and
     * non refracting atmosphere)
     * Ref: Salomons, E., Van Maercke, D., Defrance, J.,&amp;De Roo, F. (2011). The Harmonoise sound propagation model.
     * Acta acustica united with acustica, 97(1), 62-74 (section 2.5)
     */
    @Test
    public void test_harmonoise_case01() throws IOException {

        double[] attenuation = computeHarmonoiseAttenuation("case_1_4");

        //Assertion
        assertEquals(0, attenuation[0]);

    }
}
