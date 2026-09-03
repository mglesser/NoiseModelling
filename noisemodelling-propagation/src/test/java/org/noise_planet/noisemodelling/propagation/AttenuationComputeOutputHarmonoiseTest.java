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
     * Acta acustica united with acustica, 97(1), 62-74 (section 3)
     */
    @Test
    public void test_harmonoise_case01() throws IOException {

        double[] attenuation = computeHarmonoiseAttenuation("case_05");

        //Assertion
        double[] referenceExcessAttenuation = {
                6.036446469248293	,
                6.10478359908884	,
                6.10478359908884	,
                6.10478359908884	,
                6.10478359908884	,
                6.10478359908884	,
                6.10478359908884	,
                6.036446469248293	,
                5.9681093394077465	,
                5.8997722095672	,
                5.831435079726653	,
                5.694760820045559	,
                5.4897494305239185	,
                5.148063781321184	,
                4.533029612756266	,
                3.4396355353075165	,
                1.662870159453302	,
                -1.8906605922551236	,
                -10.979498861047835	,
                -4.760820045558086	,
                2.619589977220958	,
                5.626423690205012	,
                4.3963553530751724	,
                -3.8724373576309787	,
                4.191343963553532	,
                1.457858769931665	,
                4.259681093394079
        }; // plotdigitized from publication, 1/3 oct band 25Hz - 10kHz
        assertEquals(0, attenuation[0]);

    }
}
