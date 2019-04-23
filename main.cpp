#include <iostream>

#include <stdio.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <fcntl.h>


#include <cstring>
#include <fstream>

#include <assert.h>
#include <string.h>
#include <map>
#include <ctype.h>
#include <vector>
#include <cstdio>

#include <cstring>
#include <stdio.h>
#include <cmath>
#include <zconf.h>
#include<map>

#include<iostream>
using namespace std;
#include <fstream>
#include <cmath>
#include <complex>
#include <stdlib.h>
#include <stdio.h>
#define PI 3.141592654
#include<random>
#include<time.h>
#include <Eigen/Dense>
#include"fft.hpp"
#include "Fft1.h"
using namespace Eigen;
using namespace Eigen::internal;
using namespace Eigen::Architecture;
using namespace Eigen ;



void griffin_lim(MatrixXf *S);
VectorXf istft(const MatrixXcf *S,const VectorXf*win);
MatrixXcf stft(VectorXf *y,VectorXf*win);
MatrixXf compute_angle(MatrixXcf *);
VectorXf ifft(VectorXcf *V,int length);
VectorXcf fft(VectorXf *V,int length);
VectorXcf fft2(VectorXf *V,int length);
float *inv_mel_basis = new float[513*80];
float *ifft_window = new float[1024];
VectorXf ifft_wind(1024);
complex<float> complex3j(0,3);

VectorXf window_sumsquare(VectorXf *win,int n_frames);
MatrixXf compute_angle(MatrixXcf *S){
    MatrixXf out(513,316);
    //std::cout<<(*S);
    MatrixXcf S2 = *S;
    //std::cout<<complex3j.imag();
    for (int i = 0; i < S->size(); i++) {
        complex<float> tmp = *(S->data()+i);
        *(out.data()+i) = atan2(tmp.imag(),tmp.real())*1.0; //57.295779513;
    }
    //std::cout<<out;
    return out;
}

void denormalize(float *mel,int length){
    for(int i=0;i<length;i++){
        if(mel[i]>4.0){
            mel[i] = 4.0;
        }
        if(mel[i]< -4.0){
            mel[i] = -4.0;
        }
        mel[i] = 12.5*(mel[i]+4.0) - 100.0;
        mel[i] = pow(10.0,mel[i]*0.05+1.);


    }

}

MatrixXcf stft(VectorXf *y,VectorXf*win){
    MatrixXcf out(513,316);
    int win_length = 1024;
    int hop_length = 256;
    //float *y_pad = new float(81644);
    VectorXf y_pad(81664);
    y_pad.setZero();
    //std::cout<<"y:"<<std::endl;
    for (int i = 0; i < 512; i++) {
        *(y_pad.data()+i) = y->data()[512-i];
    }

    for (int i = 0; i < 80620; i++) {
        *(y_pad.data()+512+i)=y->data()[i];
    }
    for (int i = 0; i < 512; i++) {
        *(y_pad.data()+81644-i) = (*y)[80620-512+i];
    }
    MatrixXf tmp(1024,316);

    for (int i = 0; i < 316; i++) {
        tmp.col(i) = y_pad.segment(i*hop_length,1024);
    }
    MatrixXcf y_fft(1024,316);

    VectorXf V(1024);

    for (int i = 0; i < 316; i++) {
        V = (*win).cwiseProduct(tmp.col(i));
        y_fft.col(i)=fft(&V,1024);

    }
    return y_fft.topRows(513);

}




VectorXf istft(MatrixXcf *S, VectorXf *win){
    int win_length = 1024;
    int hop_length = 256;
    VectorXf out(81664);

    //out.setZero();
    //float *win = new float(win_length);
    int n_frames = S->cols();
    //std::cout<<std::endl<<"n_frames:"<<n_frames<<std::endl;
    //std::cout<<std::endl<<"S:"<<*win<<std::endl;
    int wav_len = 81644;
    VectorXf win_squa =window_sumsquare(win,316);

    //std::cout<<"win_squa"<<win_squa<<std::endl;
    for(int i=0;i<316;i++){
        int start_sam = i*hop_length;
        VectorXcf spec1 =S->col(i);
        VectorXcf spec2(spec1.size()-2);

        for(int j=0;j<spec2.size();j++){
            *(spec2.data()+j) = *(spec1.data()+spec1.size()-1-1-j);
        }
        VectorXcf spec(spec1.size()+spec2.size());
        spec<<spec1,spec2;
        //std::cout<<std::endl<<"spec1:"<<spec<<spec.size()<<std::endl;
        //VectorXf spec_iff(1024);
        VectorXf spec_iff = ifft(&spec,1024);
        if (i == 0) {
        }

        VectorXf y_tmp =  (*win).cwiseProduct(spec_iff);


        for(int j=0;j<1024;j++){

            *(out.data()+start_sam+j) = *(out.data()+start_sam+j)+*(y_tmp.data()+j);
        }


    }

    for (int i = 0; i < win_squa.size(); i++)
        if( win_squa[i]!=0)
            *(out.data()+i) = *(out.data()+i) / *(win_squa.data() +i);
    VectorXf out_new = out.segment(512,80640);

    return out_new;


}
VectorXf ifft(VectorXcf *V,int length){
    VectorXf real=V->real();
    VectorXf imag=V->imag();

    IFFT(real.data(),imag.data(),length);


    return  real;
}
VectorXcf fft(VectorXf *V,int length){
    //VectorXcf out(length);
    complex<float> complex1j(0,1);
    VectorXf imag(length);
    imag.setZero();
    VectorXf real = *V;
    FFT(real.data(),imag.data(),length);

    VectorXcf out=real+complex1j*imag;
    return  out;

}
VectorXcf fft2(VectorXf *V,int length){
    complex<float> complex1j(0,1);
    VectorXf imag(length);
    imag.setZero();
    VectorXf real = *V;


    float *vec = V->data();
    for (int i = 0; i < 16; i++) {
        //std::cout<<vec[i]<<endl;
    }
    int len = length;

    Complex *inVec = new Complex[len];
    Complex *outVec = new Complex[len];
    Complex *invert = new Complex[len];

    memset(inVec, 0, len*sizeof(Complex));
    for (int i = 0; i < len; i++)
        inVec[i].rl = vec[i];

    CFft1 t;

    t.fft(inVec, len, outVec);

    for (int i = 0; i < len; i++) {
        std::cout<<outVec[i].rl<<"  "<< outVec[i].im<<endl;
    }



    VectorXcf out=real+complex1j*imag;
    return  out;
}

VectorXf window_sumsquare(VectorXf *win,int n_frames){
    int length = 81664;
    int hop_size = 256;
    VectorXf out(length);
    VectorXf win_square(1024);
    out.setZero();
    win_square.setZero();
    for (int i = 0; i < 1024; i++) {
        *(win_square.data()+i)=*(win->data()+i) * *(win->data()+i);
    }

    for (int i = 0; i < n_frames; i++) {
        int start_i = i*hop_size;
        for(int j=0;j<1024;j++){
            if(start_i+j>length-1)
                break;
            *(out.data()+start_i+j) += *(win_square.data()+j);
        }

    }
    return out;
}

float gaussrand()
{
    static float U, V;
    static int phase = 0;
    float z;

    if(phase == 0)
    {
        U = rand() / (RAND_MAX + 1.0);
        V = rand() / (RAND_MAX + 1.0);
        z = sqrt(-2.0 * log(U))* sin(2.0 * PI * V);
    }
    else
    {
        z = sqrt(-2.0 * log(U)) * cos(2.0 * PI * V);
    }

    phase = 1 - phase;
    return z;
}


class WavWriter {
private:
    std::ofstream stream;
    template <typename T>
    void write(std::ofstream& stream, const T& t) {
        stream.write((const char*)&t, sizeof(T));
    }

    void writeWAVData(size_t buf, size_t samples, int sampleRate, int16_t channels)
    {
        int bytes=samples*channels*buf;

        stream.write("RIFF", 4);
        write<int32_t>(stream, 36 + bytes);
        stream.write("WAVE", 4);
        stream.write("fmt ", 4);
        write<int32_t>(stream, 16);
        write<int16_t>(stream, 1);                         // Format (1 = PCM)
        write<int16_t>(stream, channels);                  // Channels
        write<int32_t>(stream, sampleRate);                  // Sample Rate
        write<int32_t>(stream, sampleRate * channels * buf); // Byterate
        write<int16_t>(stream, channels * buf);            // Frame size
        write<int16_t>(stream, 8 * buf);                   // Bits per sample
        stream.write("data", 4);
        stream.write((const char*)&bytes, 4);
    }
public:
    WavWriter (const char* outFile, size_t buf,
               size_t samples, int sampleRate, int16_t channels) :
            stream (outFile, std::ios::binary) {
        writeWAVData (buf, samples, sampleRate, channels);
    }

    template <typename SampleType> void writeFrame (SampleType *buf, size_t bufsize)
    {
        stream.write ((const char*)buf, bufsize);
    }

};
int vectorToWav (VectorXf* v_in, const std::string& file, int sr)
{
    float max;

    //int sr1 = 16000;

    VectorXf v = *v_in;

    max=v.array().abs().maxCoeff();
    v=v*(1/max)*32767;

    WavWriter wv (file.c_str(), sizeof(int16_t), v.size(), sr, 1);

    for (int i=0; i<v.size(); ++i) {
        int16_t value=(int16_t)v(i);
        wv.writeFrame (&value, sizeof(int16_t));
    }

    return 0;
}



int main()
{

    const int mel_len = 316*80;
    float *mels = new float[mel_len];

    std::ifstream ifs4("mels", std::ios::binary | std::ios::in);
    if(!ifs4){
        cerr << "open file error " << endl;
        exit(0);
    }
    ifs4.read((char*)mels, sizeof(float)*mel_len);
    ifs4.close();
    std::ifstream ifs5("inv_mel_basis", std::ios::binary | std::ios::in);
    if(!ifs5){
        cerr << "open file error " << endl;
        exit(0);
    }

    ifs5.read((char*)inv_mel_basis, sizeof(float)*513*80);
    ifs5.close();

    std::ifstream ifs6("ifft_window", std::ios::binary | std::ios::in);

    ifs6.read((char*)ifft_window, sizeof(float)*1024);
    ifs6.close();
    VectorXf wind(1024);
    for (int i = 0; i < 1024; i++) {
        *(wind.data() + i) = ifft_window[i];
    }



    //denormalize(mels,mel_len);


    Eigen::Map<Matrix<float ,80,316> > mat_mel(mels);
    Eigen::Map<Matrix<float ,513,80> > mat_invmel(inv_mel_basis);



    MatrixXf linear = mat_invmel*mat_mel;
    float *p = linear.data();
    for (int i = 0; i < 513*316; i++)
        if(p[i]<1e-10)
            p[i] = 1e-10;
    //std::cout << std::endl << linear << std::endl;
    for (int i = 0; i < 513*316; i++)
        p[i] = pow(p[i] , 1.55);
    //std::cout << std::endl << linear << std::endl;
    complex<float > complex2j(0,2),complex1j(0,1);


    //std::cout << std::endl << complex2j << std::endl;

    MatrixXf rand_mat(513, 316);
    //MatrixXcf rand_matc(513, 316);
    //std::cout << std::endl << rand_matc(2,2)<<" "<<rand_matc(2,5) << std::endl;
    std::default_random_engine random(time(NULL));
    std::uniform_real_distribution<float> dis2(0.0, 1.0);
    for (int i = 0; i < rand_mat.size(); i++) {
        *(rand_mat.data() + i) = dis2(random);
    }

    rand_mat.setOnes();
    //std::cout << std::endl << rand_mat << std::endl;
    float pi = 3.141592653589793;
    MatrixXcf angles = complex2j * pi * rand_mat;
    for(int i=0;i<angles.size();i++)
        *(angles.data()+i) = exp(*(angles.data()+i));


    complex<float> complex1(1,0);
    MatrixXcf S_complex = complex1*linear;
    MatrixXcf tmp = S_complex.cwiseProduct(angles);
    VectorXf y =istft(&tmp,&wind);
    for(int i=0;i<60;i++){
        MatrixXcf y_sf   = stft(&y,&wind);

        if (i == 0) {
            //std::cout << std::endl <<  std::endl << "y_sf:"<<y_sf.col(270) << std::endl;
        }


        MatrixXf y_angle = compute_angle(&y_sf);

        angles           = complex1j * y_angle;


        for(int j=0;j<angles.size();j++)
            *(angles.data()+j) = exp(*(angles.data()+j));
        tmp = S_complex.cwiseProduct(angles);
        y = istft(&tmp,&wind);

    }
    vectorToWav (&y, "output2.wav", 16000);
/*



    //std::cout << std::endl <<  std::endl << "y" << std::endl;



    //std::cout << std::endl << angles.size() << std::endl << linear.size() << std::endl<< tmp.size() << std::endl;





    //MatrixXf tmp = linear.transpose();
    //std::cout << std::endl << angles << std::endl << angles.size() << std::endl;

    //std::cout << std::endl <<  std::endl << angles.size() << std::endl;

    //std::cout << std::endl << tmp << std::endl;
    //for (int i = 0; i < 100; i++)
        //std::cout << *(tmp.data() + i) << " ";
    //std::cout << std::endl << std::endl;


    //for (int i = 0; i < 100; i++)
        //std::cout << *(mat_mel.data() + i) << " ";
    //std::cout << std::endl << std::endl;
    /*
    for (int i = 0; i < 100; i++)
        std::cout << *(mat_invmel.data() + i) << " ";
    std::cout << std::endl << std::endl;



    for(int i=0;i<513*80;i++){
        std::cout<<inv_mel_basis[i]<<' ';
        if(i%160 == 0){
            std::cout<<endl<<i/160<<endl;
        }
    }


*/




    return 0;
}
