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
import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.noise_planet.noisemodelling.pathfinder.profilebuilder.CutProfile;

import java.util.List;

import static java.lang.Math.*;

/**
 * 2D ground profile between a source and a receiver.
 * @author Martin Glesser
 */

public class HarmonoiseGroundProfile {
    LineString profile; // 2D coordinates of ground profile vertices
    Coordinate source; // 3D coordinates of the source
    Coordinate receiver; // 3D coordinates of the receiver
    Coordinate[] vertices;

    /**
     * Initialize HarmonoiseGroundProfile object from CutProfile object.
     *
     * @param cutProfile 3D profile from source to receiver
     */
    public HarmonoiseGroundProfile(CutProfile cutProfile){
        source = cutProfile.getSource().getCoordinate();
        receiver = cutProfile.getReceiver().getCoordinate();
        // Get the whole 2D profile including ground points
        List<Integer> hullIndices = cutProfile.getConvexHullIndices(cutProfile.computePts2D());
        vertices = cutProfile.computePts2DGround(hullIndices).toArray(new Coordinate[0]);
        GeometryFactory geometryFactory = new GeometryFactory();
        profile = geometryFactory.createLineString(vertices);
    }

    /**
     * Generate a curved profile from a coordinate list, two endpoints (source and receiver) and a curvature radius.
     * Ref: Salomons, E., Van Maercke, D., Defrance, J.,&amp;De Roo, F. (2011). The Harmonoise sound propagation model.
     * Acta acustica united with acustica, 97(1), 62-74 (section 2.5)
     * Note: This implementation yield similar results to the one from CurvedProfileGenerator.applyTransformation.
     * However, it works only on the whole ground profile (from zGroundSource to zGroundReceiver).
     */
    public void computeCurvedProfile(double radius){
        // Segment Profile (second paragraph of section 2.5)
        double dsr = source.distance(receiver);
        double maxSegmentLength = min( dsr/3 , max(50, dsr/20));
        profile = (LineString) Densifier.densify(profile, maxSegmentLength);

        // Ground curvature
        double hSource = source.z;
        double hReceiver = receiver.z;
        double hm = (hSource + hReceiver) / 2;
        double c0 = 2* (hm + radius); // Eq. 77
        Complex c = new Complex(0, c0); // Eq. 76
        double xc = 0.5 * (profile.getStartPoint().getX() + profile.getEndPoint().getX());
        double yc = 0.5 * (profile.getStartPoint().getY() + profile.getEndPoint().getY()) + hm;
        Complex w0 = new Complex(xc, yc); // Eq. 75
        double deltaY = 0;
        vertices = new Coordinate[profile.getNumPoints()];
        for (int i = 0; i < profile.getNumPoints(); i++) {
            Complex w = new Complex(profile.getCoordinateN(i).getX(), profile.getCoordinateN(i).getY());
            Complex wPrim = c.multiply(w.subtract(w0)).divide(c.add(w.subtract(w0))); // Eq. 74

            // Create new coordinate with transformed z (incl. profile translation)
            if (i == 0) {
                deltaY = profile.getCoordinateN(i).getY() - wPrim.getImaginary();
                vertices[i] =
                        new Coordinate(wPrim.getReal() + xc, profile.getCoordinateN(i).getY() , profile.getCoordinateN(i).getZ());
            } else {
                vertices[i] =
                        new Coordinate(wPrim.getReal() + xc, wPrim.getImaginary() + deltaY, profile.getCoordinateN(i).getZ());
            }
        }
    }

    public Coordinate getSource() {
        return source;
    }

    public Coordinate getReceiver() {
        return receiver;
    }
}
