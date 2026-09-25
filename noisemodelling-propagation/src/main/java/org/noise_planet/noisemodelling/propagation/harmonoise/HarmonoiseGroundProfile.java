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
import org.locationtech.jts.algorithm.Angle;
import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.math.Vector2D;
import org.noise_planet.noisemodelling.pathfinder.profilebuilder.CutProfile;
import org.noise_planet.noisemodelling.pathfinder.profilebuilder.SurfaceAbsorption;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static java.lang.Math.*;

/**
 * 2D ground profile between a source and a receiver.
 * @author Martin Glesser
 */

public class HarmonoiseGroundProfile {
    private Coordinate[] vertices;
    private double[] flowResistivity;
    private Complex[][] reflectionCoefficient;

    /**
     * Initialize HarmonoiseGroundProfile object from an array of vertices.
     * (for testing purpose)
     *
     * @param vertices 2D ground profile vertices
     */
    public HarmonoiseGroundProfile(Coordinate[] vertices){
        this.vertices = vertices;
    }

    /**
     * Initialize HarmonoiseGroundProfile object from CutProfile object.
     *
     * @param cutProfile 3D profile from source to receiver
     * @param nFreq size of the frequency axis
     */
    public HarmonoiseGroundProfile(CutProfile cutProfile, int nFreq){
        this(cutProfile, nFreq, 0);
    }

    /**
     * Initialize HarmonoiseGroundProfile object from CutProfile object.
     *
     * @param cutProfile 3D profile from source to receiver
     * @param nFreq size of the frequency axis
     * @param radius curvature radius of the equivalent ground profile
     */
    public HarmonoiseGroundProfile(CutProfile cutProfile, int nFreq, double radius){
        // Get the whole 2D profile including ground points
        extractVertices(cutProfile);
        // TODO : Manage geometry densification
        // Densify per segment to update flowResistivity array
        // TODO : Manage ground curvature
        reflectionCoefficient = new Complex[getNVertices()-1][nFreq];
    }

    private void extractVertices(CutProfile cutProfile){
        List<Integer> indices = new ArrayList<>(0);
        vertices = cutProfile.computePts2DGround(0, indices).toArray(new Coordinate[0]);
        flowResistivity = indices.stream()
                .mapToDouble(i -> cutProfile.cutPoints.get(i).groundCoefficient)
                .toArray();
        // TODO : Manage flow resistivity for all cutPoints type
    }

    /**
     * Generate a curved profile from a coordinate list, two endpoints (source and receiver) and a curvature radius.
     * Ref: Salomons, E., Van Maercke, D., Defrance, J.,&amp;De Roo, F. (2011). The Harmonoise sound propagation model.
     * Acta acustica united with acustica, 97(1), 62-74 (section 2.5)
     * Note: This implementation yield similar results to the one from CurvedProfileGenerator.applyTransformation.
     * However, it works only on the whole ground profile (from zGroundSource to zGroundReceiver).
     */
    private Coordinate[] computeCurvedProfile(Coordinate[] groundProfile, double radius){
        GeometryFactory geometryFactory = new GeometryFactory();
        LineString profile = geometryFactory.createLineString(groundProfile);
        Coordinate source = groundProfile[0];
        Coordinate receiver = groundProfile[groundProfile.length-1];
        Coordinate[] vertices;
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
        return vertices;
    }

    public Coordinate[] getVertices() {
        return vertices;
    }

    public int getNVertices() {
        return vertices.length;
    }

    /**
     * @return ground profile vertices with first and last ground vertices replaced
     * respectively by the source and the receiver vertices
     */
    public Coordinate[] getVertices(int iStart, int iEnd) {
        return Arrays.copyOfRange(vertices, iStart, iEnd+1);
    }

    /**
     * Getter for vertex coordinate.
     *
     * @param i index of the vertex
     * @return coordinate of the vertex
     */
    public Coordinate getVertex(int i){
        return vertices[i];
    }

    /**
     * Return the coordinates of the image of a vertex with respect to a ground segment plane
     *
     * @param iVertex index of the vertex
     * @param iSeg index of the first point of the ground segment
     * @return image of the vertex
     */
    public Coordinate getImageVertex(int iVertex, int iSeg) {
        double x0 = vertices[iVertex].getX();
        double y0 = vertices[iVertex].getY();
        double x1 = vertices[iSeg].getX();
        double x2 = vertices[iSeg+1].getX();
        double y1 = vertices[iSeg].getY();
        double y2 = vertices[iSeg+1].getY();
        // Segment line equation: ax + by + c = 0
        double a = 1 / (x2 - x1);
        double b = 1 / (y1 - y2);
        double c = y1 / (y2 - y1) - x1 / (x2 - x1);
        // Image vertex coordinates
        double xi = x0 - 2*a * (a*x0 + b*y0 + c) / (a*a + b*b);
        double yi = y0 - 2*b * (a*x0 + b*y0 + c) / (a*a + b*b);
        return new Coordinate(xi, yi);
    }

    /**
     * Return the abscissa of a point in a local coordinate system with the origin at the
     * normal projection of the source on the ground segment plane and the abscissa along the plane
     *
     * @param iPoint index of the point
     * @param iSeg index of the ground segment
     * @param iSource index of the (secondary) source
     * @return local abscissa (or local distance)
     */
    public double getLocalAbscissa(int iPoint, int iSeg, int iSource){
        double thetaPoint;
        if (iPoint == iSeg+1) {
            thetaPoint = 0;
        } else {
            thetaPoint = Angle.angleBetween(vertices[iPoint], vertices[iSeg + 1], vertices[iSeg]);
        }
        double thetaSource = Angle.angleBetween(vertices[iSource], vertices[iSeg+1], vertices[iSeg]);
        return vertices[iSeg+1].distance(vertices[iSource]) * Math.cos(thetaSource)
                - vertices[iSeg+1].distance(vertices[iPoint]) * Math.cos(thetaPoint);
    }

    /**
     * Return the ordinate of a point in a local coordinate system with the origin at the
     * normal projection of the source on the ground segment plane and the abscissa along the plane
     * @param iPoint index of the point
     * @param iSeg index of the ground segment
     * @return local abscissa (or local height)
     */
    public double getLocalOrdinate(int iPoint, int iSeg) {
        double theta = Angle.angleBetweenOriented(vertices[iPoint], vertices[iSeg + 1], vertices[iSeg]);
        return vertices[iSeg + 1].distance(vertices[iPoint]) * Math.sin(theta);
    }

    /**
     * Return the reflexion angle between source and receiver with respect to the normal on the ground
     * segment surface
     *
     * @param iSeg iSeg index of the ground segment
     * @param iSource index of the (secondary) source
     * @param iReceiver index of the (secondary) receiver
     * @return reflexion angle in radian
     */
    public double getReflexionAngle(int iSeg, int iSource, int iReceiver) throws IllegalArgumentException {
        if (isConvexSegment(iSeg, iSource, iReceiver)) {
            throw new IllegalArgumentException("Reflexion angle cannot be computed for convex segments");
        } else {
            Vector2D reflexionPlane = new Vector2D(vertices[iSeg], vertices[iSeg + 1]);
            Vector2D imageSourceReceiverPlane = new Vector2D(this.getImageVertex(iSource, iSeg), vertices[iReceiver]);
            return Math.PI / 2 - reflexionPlane.angle(imageSourceReceiverPlane);
        }
    }

    /**
     * Check if a segment is convex according to Harmonoise definition
     *
     * @param iSeg index of the ground segment
     * @param iSource index of the (secondary) source
     * @param iReceiver index of the (secondary) receiver
     * @return true if the segment is convex
     */
    public boolean isConvexSegment(int iSeg, int iSource, int iReceiver){
        double localSourceHeight = getLocalOrdinate(iSource, iSeg);
        double localReceiverHeight = getLocalOrdinate(iReceiver, iSeg);
        return localSourceHeight < 0 || localReceiverHeight < 0;
    }

    /**
     * Get ground segment normalized ground impedance
     *
     * @param iSeg index of the ground segment
     * @param frequency frequency axis
     * @return ground impedance
     */
    public Complex[] getGroundImpedance(int iSeg, List<Double> frequency) {
        return (Complex[]) frequency.stream()
                .map(f -> SurfaceAbsorption.computeSurfaceImpedance(flowResistivity[iSeg], f))
                .toArray();
    }

    /**
     * Set the frequency dependant reflectionCoefficient for segment iSeg
     *
     * @param iSeg index of the ground segment
     * @param reflectionCoefficient spherical-wave reflection coefficient
     */
    public void setReflectionCoefficient(int iSeg, Complex[] reflectionCoefficient) {
        this.reflectionCoefficient[iSeg] = reflectionCoefficient;
    }

    /**
     * Return the reflection coefficient at a given frequency
     * @param iSeg index of the ground segment
     * @param iFreq index of the frequency
     * @return reflection coefficient
     */
    public Complex getReflectionCoefficient(int iSeg, int iFreq) {
        return reflectionCoefficient[iSeg][iFreq];
    }
}
