package org.noise_planet.noisemodelling.propagation;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.noise_planet.noisemodelling.propagation.harmonoise.HarmonoiseGroundProfile;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HarmonoiseGroundProfileTest {
    private static HarmonoiseGroundProfile profile;
    private static Coordinate[] projections;

    @BeforeAll
    public static void setUp(){
        Coordinate[] vertices = new Coordinate[] {
                new Coordinate(0,3),
                new Coordinate(5, -2),
                new Coordinate(6,1), // seg
                new Coordinate(8,2),
                new Coordinate(11,1),
                new Coordinate(11, 6)
        };
        profile = new HarmonoiseGroundProfile(vertices);
        projections = new Coordinate[] {
                new Coordinate(2, -1),
                new Coordinate(4, 0),
                new Coordinate(6,1),
                new Coordinate(8,2),
                new Coordinate(10,3),
                new Coordinate(12,4)
        };
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
        int i = 4;get
        d = profile.getLocalAbscissa(i,2, 2);
        assertEquals(profile.getVertex(2).distance(projections[i]), d);
    }
}
