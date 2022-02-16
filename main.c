#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "tf1d.c"
#include "cpgplot.h"

/* 'M_PI' is most certainly defined in 'math.h' on linux but is not
   actually part of the C standard library. */

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

FILE* file;

/* ================================ 1/
/* Defining the function for the FT */
/* ================================ */

float func(float t, float nu, int cos_flag)
    {
        float y;
        if(cos_flag == 1) y = cos(2*M_PI*nu*t);
        if(cos_flag == 2) y = cos(2*M_PI*nu*t*t);
        if(cos_flag == 3) y = cos(2*M_PI*nu*t*t*t);
        else y = cos(2*M_PI*nu*powf(t,cos_flag)); /* for powers above 3 we use the not so efficient pow */
        return y;
    }

float gauss(float t, float sigma, float mu)
    {
        float y;
        y = (1 / sqrt(2*M_PI*sigma*sigma)) * expf(-((t-mu)*(t-mu))/(2*sigma*sigma));
        return y;
    }

int main()
{
    /* ===================== */
    /* Variable declarations */
    /* ===================== */

    int freq_flag;
    int trunc_flag;
    int read_flag;
    int cos_flag;
    float h; /* resolution step */
    float dt; /* truncation interval */
    int N;

    file = fopen("conf.txt","r");
    fscanf(file,"%d",&freq_flag);
    fscanf(file,"%d",&trunc_flag);
    fscanf(file,"%d",&read_flag);
    fscanf(file,"%d",&cos_flag);
    fscanf(file,"%f",&h);
    fclose(file);

    printf("%d\t%d\t%d\t%d\t%f\n",freq_flag,trunc_flag,read_flag,cos_flag,h);

    if(read_flag == 0) N = 1024;
    else N = 32768;

    int n_line = 0;
    float tmax = (N-1)*h;
    float nu1 = 1.0 / (((1.0*N)/20)*h); /* in Hz */;
    float nu2 = 4.0 / (((1.0*N)/20)*h); /* in Hz */;
    float tp = 0.; /* truncation point */
    float gmin = 0, gmax = 0;
    float Wmin = 0, Wmax = 0;
    float *t; /* time domain */
    float *im;
    float *gaussian;
    float *freq; /* frequency domain */
    float **W = malloc(N*sizeof(float*));
    for (unsigned int i = 0; i < N; i++) W[i] = malloc((N/2)*sizeof(float));
    float **g = malloc(N*sizeof(float*));
    for (unsigned int i = 0; i < N; i++) g[i] = malloc(N*sizeof(float));
    float buffer1, buffer2; /* buffers for file reading */

    t = (float *)malloc(N*sizeof(float));
    im = (float *)malloc(N*sizeof(float));
    gaussian = (float *)malloc(N*sizeof(float));
    freq = (float *)malloc((N/2)*sizeof(float));


    float *tab;
    float xmin,xmax;
    float ymin,ymax;
    float zmin,zmax;
    float tr[6];

    tab = (float *)malloc(N*(N/2)*sizeof(float));

    file = fopen("piano2.wav.txt","r");

    while(fscanf(file,"%f %f\n",&buffer1,&buffer2) != EOF) 
    {
        n_line++; /* we read the file once to find the amount of line */
    }

    fclose(file);

    float **full_file = malloc(n_line*sizeof(float*));
    for (unsigned int i = 0; i < n_line; i++) full_file[i] = malloc(2*sizeof(float));
    
    for(unsigned int ix = 0; ix<n_line; ix++) 
    {
        full_file[ix][0] = 0.;
        full_file[ix][1] = 0.;
    }

    n_line = 0;

    file = fopen("piano2.wav.txt","r");
    
    while(!feof(file))
    {
        fscanf(file,"%f %f\n",&full_file[n_line][0],&full_file[n_line][1]); /* we read the file a second time to populate the matrix */
        n_line++;
    }

    fclose(file);

    /* ================================= */
    /* Time and frequency Initialisation */
    /* ================================= */

    for(unsigned int ix = 0; ix < N; ix++)
    {
        if(read_flag == 0) t[ix] = ix*h;
        else t[ix] = full_file[ix][0];
    }

    for(unsigned int ix = 0; ix < N/2; ix++)
        {
            freq[ix] = 1.*ix / (N * h);
        }

    for(unsigned int tx = 0; tx < N; tx++)
    {
        for(unsigned int ix = 0; ix < N; ix++)
        {
            g[tx][ix] = 0.;
        }
        for(unsigned int ix = 0; ix < N/2; ix++)
        {
            W[tx][ix] = 0.;
        }
    }

    dt = (t[N-1]/20.0); /*the truncation interval is defined dynamically */

    if(N <= 4096) cpgbeg(0,"/xw",1,3);
    else cpgbeg(0,"/xw",1,2);

    /* ========= */
    /* Main Loop */
    /* ========= */
    for(unsigned int tx = 0; tx < N; tx ++)
    {

        /* ============== */
        /* Initialisation */
        /* ============== */

        for(unsigned int ix = 0; ix < N; ix++)
        {
            if(read_flag == 0)
            {
                if(freq_flag == 1)
                {
                    if(ix <= N/2)
                    {
                        g[tx][ix] = func(t[ix], nu1, cos_flag);
                    }
                    else
                    {
                        g[tx][ix] = func(t[ix], nu2, cos_flag);
                    }
                }
                else g[tx][ix] = func(t[ix],nu1, cos_flag);
            }
            else
            {
                g[tx][ix] = full_file[ix][1];
            }
            im[ix] = 0.;
            gaussian[ix] = gauss(t[ix],dt,tp);
        }

        /* ========== */
        /* Truncation */
        /* ========== */

        tp = t[tx];

        if(trunc_flag == 0)
        {
            for(unsigned int ix = 0; ix < N; ix++)
            {
                if(tp-dt <= t[ix] && t[ix] <= tp+dt){} /* XOR gate here idk how to make one */
                else g[tx][ix] = 0.;
            }
        }
        else
        {
            for(unsigned int ix = 0; ix < N; ix++) g[tx][ix] = g[tx][ix] * gaussian[ix];
        }

        for(unsigned int ix = 0; ix < N; ix++)
        {
            if(g[tx][ix] < gmin) gmin = g[tx][ix];
            if(g[tx][ix] > gmax) gmax = g[tx][ix];
        }

        /* ====================== */
        /* Function Data graphing */
        /* ====================== */

        cpgask(0);
        cpgenv(0,t[N-1],gmin,gmax,0,1);
        cpglab("Time","Amplitude","Base signal");
        cpgline(N,t,g[tx]);

        /* == */
        /* TF */
        /* == */

        tf1d(g[tx],im,N,1);

        for(unsigned int ix = 0; ix < N/2; ix++)
        {
            W[tx][ix] = sqrt(g[tx][ix]*g[tx][ix]+im[ix]*im[ix]);
            //freq[ix] = 1.*ix / (N * h);
        }

        for(unsigned int ix = 0; ix < N/2; ix++)
        {
            if(W[tx][ix] < Wmin) Wmin = W[tx][ix];
            if(W[tx][ix] > Wmax) Wmax = W[tx][ix];
        }

        /* ============ */
        /* FFT graphing */
        /* ============ */

        cpgask(0);
        cpgenv(freq[0],freq[(N/2)-1],Wmin,Wmax,0,1);
        cpglab("Frequency","Amplitude","TF of the base signal");
        cpgline(N/2,freq,W[tx]);

        /* ============================ */
        /* Spectrogram runtime graphing */
        /* ============================ */

        //we disable the runtime spectogram for large N numbers because of O(N^3) complexity

        if(N <= 4096)
        {
            for(unsigned int tx = 0; tx < N; tx++)
            {
                for(unsigned int ix = 0; ix < N/2; ix++)
                {
                    tab[tx*(N/2)+ix] = W[tx][ix];
                }
            }
            
            /* ====== */
            /* Limits */
            /* ====== */

            xmin = freq[0]; xmax = freq[(N/2)-1];
            ymin = t[0]; ymax = t[N-1];
            zmin = Wmin; zmax = Wmax;

            /* =========== */
            /* User Matrix */
            /* =========== */

            tr[0]=xmin;
            tr[1]=(xmax-xmin)/(N/2);
            tr[2]=0;
            tr[3]=ymin;
            tr[4]=0;
            tr[5]=(ymax-ymin)/N;

            cpgask(0);
            cpgenv(xmin, xmax, ymin, ymax, 0, 1);
            cpglab("Frequency","Time","Spectrogram of the base signal");
            cpgimag(tab, (N/2), N, 1, (N/2), 1, N, zmin, zmax, tr);
            
            if(N <= 2048) usleep(10000);
        }
    }
    cpgask(0);
    cpgend();

    /* ================= */
    /* Final Spectrogram */
    /* ================= */
    
    for(unsigned int tx = 0; tx < N; tx++)
    for(unsigned int ix = 0; ix < N/2; ix++)
    {
        tab[tx*(N/2)+ix] = W[tx][ix];
    }
    
    /* ====== */
    /* Limits */
    /* ====== */

    xmin = freq[0]; xmax = freq[(N/2)-1];
    ymin = t[0]; ymax = (N-1)*h;
    zmin = Wmin; zmax = Wmax;

    /* =========== */
    /* User Matrix */
    /* =========== */

    tr[0]=xmin;
    tr[1]=(xmax-xmin)/(N/2);
    tr[2]=0;
    tr[3]=ymin;
    tr[4]=0;
    tr[5]=(ymax-ymin)/N;

    cpgbeg(0,"/xw",1,1);
    cpgenv(xmin, xmax, ymin, ymax, 0, 1);
    cpgimag(tab, (N/2), N, 1, (N/2), 1, N, zmin, zmax, tr);
    cpglab("Frequency nu (in Hz)", "Time", "Spectrogram of the base signal");

    cpgask(1);
    cpgend();
}
