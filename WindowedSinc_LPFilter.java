import ij.*;
import ij.plugin.PlugIn;
import ij.process.*;
import ij.gui.*;
import ij.text.*;
import ij.io.*;
import java.awt.*;
import java.util.*;
import java.lang.Math;

/* Implementation of windowed-sinc low-pass filter appeared in Chapter 13 of */
/* "The Scientist and Engineer's Guide to Digital Signal Processing" */
/* by Steven W. Smith, Ph.D. California Technical Publishing */
/* ISBN 0-9660176-3-3 (1997) */
/* http://www.dspguide.com/  */
/* Perform low-pass filtering against image stack regarding time axis */
/* written by k.yoshida.1 at mark bham.ac.uk */
/* 2005/10/21: First Version */
/* 2006/08/10: Second Version, high-pass filter added, 32-bit image output */
/*             Note that now the kernel size is specified by its half-width. */


public class WindowedSinc_LPFilter implements PlugIn {
    static int khw=20;
    static double fl=0.25;
    static boolean showFilter;
    static boolean smooth_ends=true;
    int width, height;
    double PI=Math.PI;
    
    public void run(String arg) {
	ImagePlus imp=WindowManager.getCurrentImage();
	if (imp==null || imp.getStackSize()==1 || imp.getType()==ImagePlus.COLOR_RGB) {
		IJ.error("This plugin requires a non-RGB stack");
		return;
	}
	GenericDialog tgd = new GenericDialog("Windowed-sinc Low-Pass filter", IJ.getInstance());
	tgd.addNumericField("kernel ±half-width:", khw, 0);
	tgd.addNumericField("Low-Pass Cutoff frequency (0-1.0)x¹:", fl, 4);
	tgd.addCheckbox("Show Filter", showFilter);
	tgd.addCheckbox("Smooth The Both Ends With Shorter Kernel", smooth_ends);
	tgd.showDialog();
	if (tgd.wasCanceled()) 
	    return;
	khw=(int)tgd.getNextNumber();
	fl=(double)tgd.getNextNumber();
	showFilter = tgd.getNextBoolean();
	smooth_ends = tgd.getNextBoolean();
	if (2*khw+1>=imp.getStackSize()) {
		IJ.error("kernel half-width must be less than half the stack size");
		return;
	}

	width=imp.getWidth();
	height=imp.getHeight();
	int nSlices = imp.getStackSize();
	ImageStack bis=imp.getStack();

	ImageProcessor bip=bis.getProcessor(1).convertToFloat();
	ImageStack ois = new ImageStack(width, height);

	for (int j=1; j<=khw; j++) {
	    IJ.showStatus("Processing "+j+"/"+nSlices+"");
	    IJ.showProgress(j, nSlices);
	    if (smooth_ends) {
		double[] eh=calcLowPassKernel(j-1,0.5*fl);
		FloatProcessor fip=convolveFrames(bis, j, j-1, eh);
		ois.addSlice("head"+j+"", fip.duplicate()); 
	    } else {
	    	ImageProcessor fip=bis.getProcessor(j).convertToFloat();
		ois.addSlice("head"+j+"", fip.duplicate()); 

	    }
	}

	double[] h=calcLowPassKernel(khw,0.5*fl);

	if (showFilter) {
	    double[] xValues=new double[2*khw+1];
	    for (int i=0;i<=2*khw;i++)
		xValues[i] = i-khw;
	    Plot plot = new Plot("Filter Kernel","Time","", xValues, h);
	    plot.draw();
	    plot.show();
	}

	for (int j=khw+1; j<=nSlices-khw; j++) {
	    IJ.showStatus("Processing "+j+"/"+nSlices+"");
	    IJ.showProgress(j, nSlices);
	    FloatProcessor fip=convolveFrames(bis, j, khw, h);
	    ois.addSlice("middle"+j+"", fip.duplicate()); 
	}

	for (int j=nSlices-khw+1; j<=nSlices; j++) {
	    IJ.showStatus("Processing "+j+"/"+nSlices+"");
	    IJ.showProgress(j, nSlices);
	    if (smooth_ends) {
		double[] eh=calcLowPassKernel(nSlices-j, 0.5*fl);
		FloatProcessor fip=convolveFrames(bis, j, nSlices-j, eh);
		ois.addSlice("tail"+j+"", fip.duplicate()); 
	    } else {
		ImageProcessor fip=bis.getProcessor(j).convertToFloat();
		ois.addSlice("foot"+j+"", fip.duplicate()); 
	    }
	}

	ImagePlus oimp=new ImagePlus(imp.getShortTitle()+"_WS"+khw+"_"+fl+"Hz.tif",ois);
	oimp.show();
    }

    double[] calcLowPassKernel(int khw, double fc) {
	double[] a=new double[2*khw+1];
	if (khw==0) {
	    a[0]=1;
	    return a;
	}

	double suma=0;
	for (int i=0;i<=2*khw;i++) {
	    if (i==khw) {
		a[i]=2*PI*fc;
	    } else {
		a[i]=Math.sin(2*PI*fc*(i-khw))/(i-khw);
	    }
	    //	    IJ.log("a[i]="+a[i]);
	    a[i]=a[i]*(0.42-0.5*Math.cos(PI*i/khw) + 0.08*Math.cos(2*PI*i/khw));
	    suma+=a[i];

	}
	for (int i=0;i<=2*khw;i++) {
	    if (suma!=0) {
		a[i]=a[i]/suma;
	    } else {
		a[i]=0;
	    }
	}
	return a;
    }
    
    FloatProcessor convolveFrames (ImageStack bis, int j, int khw, double[] h) {
	FloatProcessor fip=new FloatProcessor(width, height);
	float[] fpixels=(float[])fip.getPixels();

	for (int i=0; i<=2*khw; i++) {
	    ImageProcessor bip=bis.getProcessor(j-i+khw).convertToFloat();
	    float[] bpixels = (float[]) bip.getPixels();
	    for (int x=0;x<width;x++) {
		for (int y=0;y<height;y++) {
		    int k=width*y+x;
		    if (i==0)
			fpixels[k]=0;
		    fpixels[k]+=((double)bpixels[k])*h[i];
		}
	    }

	}
	return fip;
    }


}