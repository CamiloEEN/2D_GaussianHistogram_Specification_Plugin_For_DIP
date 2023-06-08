import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.filter.*;

public class Filter_GaussianHistTransform2D implements PlugInFilter {
    ImagePlus imp;
    
    static final int L=256;
    double a=128.0; //valor esperado
    double u=30.0; //varianza
    double k1=1.0;//Parametro 1 escala media
    double k2=1.0;//Parametro 2 escala varianza
    static final int w=1;//Determina el tamaño de la vecindad de un pixel
        
    public double DesiredGaussianDistribution(int m, int n){
        double H=1/(2*Math.PI*k2*k2*u*u);
        H=H*Math.exp( -1.0*( (m-k1*a)*(m-k1*a)+(n-k1*a)*(n-k1*a) )/(2.0*k2*k2*u*u) );
        
        return H;
    }
    
    //La siguiente función no se usa!!!!!!!!!!
    public double desiredCumProb(int m){
        
        double C=0;
        
        for(int i=0;i<m;i++)
        {
            for(int j=0;j<L;j++)
            {
                C+=DesiredGaussianDistribution(i,j);
            }
        }
        
        return C;
    }

    //Obtiene la parte entera (Floor function)
    public int intPart(double num){
        
        double fPart = num % 1;
        int iPart = (int)(num-fPart);
        
        return iPart;
        
    }
    
    public int phi(int m, int n, int pixelm, int pixeln)
    {
       
        if(m==pixelm && n==pixeln)
        {
            return 1;
        }
        
        return 0;
    }
    
    public int scrw(int [][] image, int width, int height, int m, int n){
        int sum=0;
        
        for(int i=0;i<height;i++)
        {
            for(int j=0;j<width;j++)
            {
                for(int k=-intPart(w/2.0);k<=intPart(w/2.0);k++)
                {
                    for(int l=-intPart(w/2.0);l<=intPart(w/2.0);l++)
                    {
                        if((i+k>=0 && j+l>=0) && (i+k<height && j+l<width) )
                        {
                            sum=sum+phi(m, n, image[i][j], image[i+k][j+l])*(Math.abs(m-n)+1);
                        }
                    }
                }
            }
        }
        
        return sum;
    }
    
    @Override
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        GenericDialog gd = new GenericDialog("Parametros");
                gd.addNumericField("k1: ", k1, 1);
                gd.addNumericField("k2: ", k2, 1);
                gd.showDialog();
                if (gd.wasCanceled()) 
                    return DONE;
                
                k1=gd.getNextNumber();
                k2=gd.getNextNumber();
                
		return DOES_8G;
    }

    @Override
    public void run(ImageProcessor ip) {
        
       
        int width = ip.getWidth();
        int height= ip.getHeight();
        double[][] histogram2D= new double[L][L];//intialized with zero
        int[][] pixelsMInput= new int[height][width];
        //int[][] pixelsMOutput= new int[height][width];
        byte[] pixels= (byte[]) ip.getPixels();
        double[] acumHistogram= new double[L];
        double[] desiredAcumHistogram = new double [L];
        int[] transform = new int [L];
        int cursor=0;
       
        
        //convertir imagen a matrix
        for(int y=0;y<height;y++)
                {
                    for(int x=0;x<width;x++)
                    {
                        pixelsMInput[y][x]=pixels[cursor++]&0xff;
                    }
                }
        
        //Obtención de probabilidad cumulativa de la imagen deseada
        double acumTotal=0;
        for(int i=0;i<L;i++)
        {
            for(int j=0;j<L;j++)
            {
                acumTotal+=DesiredGaussianDistribution(i,j);
            }
            desiredAcumHistogram[i] = acumTotal;
        }
        
        //Obtener histograma 2D
        int scrTotal=0;
        for(int i=0;i<L;i++)
        {
            for(int j=0;j<L;j++)
            {// Se calcula el histograma no normalizado junto al factor de normalización para ejecutar scrw
                //la menor cantidad de veces (por eficiencia)
                histogram2D[i][j]=scrw(pixelsMInput, width, height, i, j);
                scrTotal=scrTotal+(int)histogram2D[i][j];
            }
        }
        
        //Normalización del histograma 2D
        for(int i=0;i<L;i++)
        {
            for(int j=0;j<L;j++)
            {
                histogram2D[i][j]=histogram2D[i][j]/scrTotal;
            }
        }
        
        //Obtención de probabilidad cumulativa de la imagen de entrada
        acumTotal=0;
        for(int i=0;i<L;i++)
        {
            for(int j=0;j<L;j++)
            {
                acumTotal+=histogram2D[i][j];
            }
            acumHistogram[i] = acumTotal;
        }
        
        //Obtención del arreglo de la transformación
        double min=2;
        double diff=0;
        for(int m=0;m<L;m++)
        {
            min=2;
            for(int n=0;n<L;n++)
            {
                diff=Math.abs(acumHistogram[m]-desiredAcumHistogram[n]);
                if(min> diff)
                {
                    min=diff;
                    transform[m]=n;
                }
            }
        }
        
        //Transformación de la imagen de entrada a la imagen deseada
        for(int y=0;y<height;y++)
        {
            for(int x=0;x<width;x++)
            {
                pixelsMInput[y][x]=transform[(pixelsMInput[y][x])];
            }
        }
        
        
        //convertir matrix a arreglo de pixeles
        cursor=0;
        for(int y=0;y<height;y++)
                {
                    for(int x=0;x<width;x++)
                    {
                        pixels[cursor++]=(byte) pixelsMInput[y][x];
                    }
                }
        
        
        
        imp.updateAndDraw();
    }

}
