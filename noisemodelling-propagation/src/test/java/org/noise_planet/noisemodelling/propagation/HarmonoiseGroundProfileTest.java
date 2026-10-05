package org.noise_planet.noisemodelling.propagation;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.noise_planet.noisemodelling.pathfinder.profilebuilder.CutProfile;
import org.noise_planet.noisemodelling.propagation.harmonoise.HarmonoiseGroundProfile;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.noise_planet.noisemodelling.propagation.AttenuationComputeOutputHarmonoiseTest.loadCutProfile;

public class HarmonoiseGroundProfileTest {
    private static HarmonoiseGroundProfile profile;
    private static Coordinate[] projections;

    @BeforeAll
    public static void setUp() throws IOException {
        //Get test data
        URL url = HarmonoiseGroundProfileTest.class.getResource("harmonoise/testProfile.json");
        CutProfile cutProfile;
        Assertions.assertNotNull(url);
        try(InputStream inputStream = url.openStream()) {
            ObjectMapper mapper = new ObjectMapper();
            cutProfile = mapper.readValue(inputStream, CutProfile.class);
        } catch (IOException e) {
            throw new IOException(e);
        }
        // Create ground profile
        profile = new HarmonoiseGroundProfile(cutProfile, 1);
        // Projection of the vertices on the line passing through the third segment
        projections = new Coordinate[] {
                new Coordinate(2, -1),
                new Coordinate(4, 0),
                new Coordinate(6,1),
                new Coordinate(8,2),
                new Coordinate(10,3),
                new Coordinate(14,5)
        };
    }

    @Test
    public void instantiationTest(){
        // Expected profile coordinates
        Coordinate[] expectedVertices = new Coordinate[] {
                new Coordinate(0,3),
                new Coordinate(5, -2),
                new Coordinate(6,1), // seg
                new Coordinate(8,2),
                new Coordinate(11,1),
                new Coordinate(13, 7)
        };
        for (int i = 0; i < expectedVertices.length; i++) {
            assertEquals(expectedVertices[i], profile.getVertex(i));
        }
    }

    @Test
    public void getImageVertexTest(){
        // Image of vertices wrt third segment
        Coordinate[] expectedImages = new Coordinate[] {
                new Coordinate(4,-5),
                new Coordinate(3, 2),
                new Coordinate(6,1), // seg
                new Coordinate(8,2),
                new Coordinate(9,5),
                new Coordinate(15, 3)
        };
        for (int i = 0; i < expectedImages.length; i++) {
            assertEquals(expectedImages[i], profile.getImageVertex(i,2));
        }
        // Image of source wrt first segment
        assertEquals(new Coordinate(0,-7), profile.getImageVertex(0,0));
        // Image of receiver wrt last segment
        assertEquals(new Coordinate(13,-5),
                profile.getImageVertex(profile.getNVertices()-1, profile.getNVertices()-2));
    }

    @Test
    public void getLocalAbscissaTest(){
        // Input data
        Coordinate origin = new Coordinate(2,-1);
        double d;
        // Test profile
        for (int i = 0; i < profile.getNVertices(); i++) {
             d = profile.getLocalAbscissa(i,2, 0);
            assertEquals(origin.distance(projections[i]), d, 0.001);
        }
        // Test side effects - source = first point of the segment
        int i = 4;
        d = profile.getLocalAbscissa(i,2, 2);
        assertEquals(profile.getVertex(2).distance(projections[i]), d);
        // ordinate of source wrt to first segment
        assertEquals(0, profile.getLocalAbscissa(0,0, 0));
        // ordinate of receiver wrt to last segment
        assertEquals(2,
                profile.getLocalAbscissa(profile.getNVertices()-1, profile.getNVertices()-2, profile.getNVertices()-2), 0.001);
    }

    @Test
    public void getLocalOrdinateTest(){
        // Input data
        int[] signs = new int[]{1, -1, 1, 1, -1, 1};
        double h;
        // Test
        for (int i = 0; i < profile.getNVertices(); i++) {
            h = profile.getLocalOrdinate(i,2);
            assertEquals(signs[i] * Math.sqrt(
                    Math.pow(projections[i].getX()-profile.getVertex(i).getX(),2)
                            + Math.pow(projections[i].getY()-profile.getVertex(i).getY(),2)
            ), h, 0.001);
        }
        // ordinate of source wrt to first segment
        assertEquals(5, profile.getLocalOrdinate(0,0));
        // ordinate of receiver wrt to last segment
        assertEquals(6, profile.getLocalOrdinate(profile.getNVertices()-1, profile.getNVertices()-2));
    }

    @Test
    public void getReflectionAngleTest() throws IllegalArgumentException {
        // Real source to real receiver through third segment
        double angle = profile.getReflexionAngle(2, 0, profile.getNVertices()-1);
        assertEquals(1.107, angle, 0.001);
        // Convex segment
        try {
            double angle1 = profile.getReflexionAngle(2, 1, profile.getNVertices() - 1);
        } catch (IllegalArgumentException e) {
            assertEquals("Reflexion angle cannot be computed for convex segments", e.getMessage());
        }
        // Real source to real receiver through first segment
        angle = profile.getReflexionAngle(0, 0, profile.getNVertices()-1);
        assertEquals(0.748, angle, 0.001);
        // Real source to real receiver through last segment
        angle = profile.getReflexionAngle(profile.getNVertices()-2, 0, profile.getNVertices()-1);
        assertEquals(1.019, angle, 0.001);
    }

    @Test
    public void isConvexSegmentTest(){
        boolean[] expected = new boolean[] {false, true, false, true, false};
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], profile.isConvexSegment(i, 0, profile.getNVertices()-1));
        }
    }


    /**
     * Compare the Cnossos and Harmonoise ground curvature implementations.
     * Note: both implementations yield similar results. The Cnossos implementation can be used on portion of profiles
     * while the Harmonoise one shall be used only on whole profiles (from zGround_source to zGround_receiver).
     */
    @Test
    public void compareGroundCurvatureTest() throws IOException {
        // Generate profile
        CutProfile cutProfile = loadCutProfile("case_07");
        // Get curved profile through Cnossos implementation
        List<Coordinate> curvedProfile = cutProfile.computePts2D(true);
        // Get curved profile through Harmonoise implementation
        double d = cutProfile.getSource().coordinate.distance(cutProfile.getReceiver().coordinate);
        double radius = Math.max(1000, 8 * d);
        HarmonoiseGroundProfile harmonoiseProfile = new HarmonoiseGroundProfile(cutProfile, 1, radius);
    }
}
