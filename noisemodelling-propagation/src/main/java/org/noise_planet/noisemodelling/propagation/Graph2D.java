package org.noise_planet.noisemodelling.propagation;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Locale;

public class Graph2D {
    private final int dWidth; // Display width
    private final int dHeight; // Display height
    private final int margin;
    private final int textMargin;
    private final double[] x;
    private final double[][] y;
    private final String[] colors = new String[]{
            "rgb(232,126,77)", // orange
            "rgb(130,166,173)", // dark blue-green
    };
    private double xMin;
    private double xMax;
    private double yMin;
    private double yMax;
    private double[] xTicks;
    private double[] yTicks;

    /**
     * Create Graph2D object with default parameters
     *
     * @param args data to be plotted (x, y1, y2, ...)
     */
    public Graph2D(double[]... args) {
        dWidth = 400;
        dHeight = 300;
        margin = 20;
        textMargin = 20;
        this.x = args[0];
        this.y = Arrays.copyOfRange(args, 1, args.length);
        axis(50, 10000, -30, 15);
        xTicks(new double[]{100, 1000, 10000});
        yTicks(new double[]{-30, -20, -10, 0, 10});
    }

    /**
     * Plot the data in a .svg file
     *
     * @param path svg file path
     */
    public void draw(String path){
        try (FileWriter writer = new FileWriter(path)){
            // Header
            StringBuilder sb = new StringBuilder(String.format("<svg height=\"%d\" width=\"%d\">\n",
                    dHeight, dWidth));
            // White background
            sb.append(String.format(
                    Locale.ROOT,
                    "<rect width=\"%d\" height=\"%d\" x=\"0\" y=\"0\" style=\"fill:white; stroke:none\"/>\n",
                    dWidth, dHeight)
            );
            // Axis box
            sb.append(String.format(
                    Locale.ROOT,
                    "<rect width=\"%d\" height=\"%d\" x=\"%d\" y=\"%d\" style=\"fill:none; stroke-width:1;stroke:black\"/>\n",
                    dWidth -2*margin-textMargin, dHeight -2*margin-textMargin, margin+textMargin, margin)
            );
            // Ticks
            drawTicks(sb);
            plot(sb);
            legend(sb);
            // Footer
            sb.append("</svg> \n");
            // Write file
            writer.write(sb.toString());
        } catch (IOException e){
            System.out.println(
                    "An error occurred while writing"
                            + " to the file: " + e.getMessage());
        }
    }

    /**
     * Add ticks and ticks labels to the svg string
     *
     * @param sb svg string StringBuilder
     */
    private void drawTicks(StringBuilder sb){
        for (double yTick : yTicks) {
            double y = getYLin(yTick);
            sb.append(String.format(
                    Locale.ROOT,
                    "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" style=\"fill:none; stroke-width:1; stroke:black\"/>\n",
                    (double) margin + textMargin, y, (double) margin + textMargin + 5, y)
            );
            sb.append(String.format(
                    Locale.ROOT,
                    "<text x=\"%f\" y=\"%f\" text-anchor=\"end\" style=\"stroke:none; fill:black; font-size:smaller\">%.0f</text>\n",
                    (double) margin + textMargin - 5, y, yTick)
            );
        }
        for (double xTick : xTicks) {
            double x = getXLog(xTick);
            sb.append(String.format(
                    Locale.ROOT,
                    "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" style=\"fill:none; stroke-width:1; stroke:black\"/>\n",
                    x, (double) dHeight - margin - textMargin, x, (double) dHeight - margin - textMargin - 5)
            );
            sb.append(String.format(
                    Locale.ROOT,
                    "<text x=\"%f\" y=\"%f\" text-anchor=\"middle\" style=\"stroke:none; fill:black; font-size:smaller\">%.0e</text>\n",
                    x, (double) dHeight - margin - textMargin + 15, xTick)
            );
        }
    }

    /**
     * Add data series to the svg string
     *
     * @param sb svg string StringBuilder
     */
    private void plot(StringBuilder sb){
        for (int j = 0; j < y.length; j++) {
            for (int i = 1; i < x.length; i++) {
                sb.append(String.format(
                        Locale.ROOT,
                        "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" style=\"fill:none; stroke-width:1; stroke:%s\"/>\n",
                        getXLog(x[i - 1]), getYLin(y[j][i - 1]), getXLog(x[i]), getYLin(y[j][i]), colors[j])
                );
            }
        }
    }

    private void legend(StringBuilder sb){
        sb.append(String.format(
                Locale.ROOT,
                "<text x=\"%f\" y=\"%d\" text-anchor=\"middle\" style=\"stroke:none; fill:black\">Frequency [Hz]</text>\n",
                margin + textMargin + (dWidth -2*margin-textMargin)/2.0, dHeight - 5)
        );
        sb.append(String.format(
                Locale.ROOT,
                "<text x=\"%d\" y=\"%f\" transform=\"rotate(-90, %d, %f)\" text-anchor=\"middle\" style=\"stroke:none; fill:black\">Excess attenuation [dB]</text>\n",
                5, margin + (dHeight -2*margin-textMargin)/2.0, 12, margin + (dHeight -2*margin-textMargin)/2.0)
        );
    }

    /**
     * Return the y coordinate corresponding to yValue along a linear scale
     *
     * @param yValue y value of the data point
     * @return y coordinate in the graph
     */
    private double getYLin(double yValue){
        return dHeight - (margin + textMargin) - (dHeight -2*margin-textMargin) * (yValue - yMin) / (yMax - yMin);
    }

    /**
     * Return the x coordinate corresponding to xValue along a log10 scale
     *
     * @param xValue x value of the data point
     * @return x coordinate in the graph
     */
    private double getXLog(double xValue){
        return margin + textMargin + (dWidth-2*margin-textMargin) * (Math.log10(xValue) - Math.log10(xMin))
                / (Math.log10(xMax) - Math.log10(xMin));
    }

    /**
     * Set axis limits
     *
     * @param xMin x axis min value
     * @param xMax x axis max value
     * @param yMin y axis min value
     * @param yMax y axis max value
     */
    public void axis(double xMin, double xMax, double yMin, double yMax){
        this.xMin = xMin;
        this.xMax = xMax;
        this.yMin = yMin;
        this.yMax = yMax;
    }

    /**
     * Set x-axis ticks
     *
     * @param ticks x-axis ticks
     */
    public void xTicks(double[] ticks){
        xTicks = ticks;
    }

    /**
     * Set the y-axis ticks
     * @param ticks y-axis ticks
     */
    public void yTicks(double[] ticks){
        yTicks = ticks;
    }


}
