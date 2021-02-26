using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

using System.Numerics;

namespace FFT
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();

            comboBoxFFTLength.Text = "256";
        }

        struct PointD                               // структура, для хранения точек сигнала типв double
        {
            public double X;
            public double Y;
        };

        void DrawSignal(PointD[] signal)            // функция рисования точек сигнала
        {
            double[] x = new double[signal.Length];
            double[] y = new double[signal.Length];

            for (int i = 0; i < x.Length; i++)
            {
                x[i] = signal[i].X;
                y[i] = signal[i].Y;
            }

            chart1.Series[0].Points.DataBindXY(x, y);
        }

        void DrawSpectrum(double[] spectrum)        // функция рисования точек спектра
        {
            chart2.Series[0].Points.DataBindXY(Get_t(spectrum), spectrum);
        }

        void DrawSpectrumClear(double[] spectrum)   // функция рисования точек очищенного спектра
        {
            chart2.Series[1].Points.DataBindXY(Get_t(spectrum), spectrum);
        }

        void DrawSignalResult(double[] signal)      // функция рисования точек восстановленного сигнала
        {
            chart3.Series[0].Points.DataBindXY(Get_t(signal), signal);
        }

        void DrawSignalOriginal(PointD[] signal)    // функция рисования точек исходного сигнала
        {
            double[] x = new double[signal.Length];
            double[] y = new double[signal.Length];

            for (int i = 0; i < x.Length; i++)
            {
                x[i] = Convert.ToDouble(i); 
                y[i] = signal[i].Y;
            }

            chart3.Series[1].Points.DataBindXY(x, y);
        }

        double[] Get_t(double[] arr)                // функция генерации точек t
        {
            double[] t = new double[arr.Length];
            for (int i = 0; i < arr.Length; i++)
            {
                t[i] = Convert.ToDouble(i);
            }

            return t;
        }

        PointD[] GeneratingPoints(double A, double f, double Phi, int length, bool b)  // функуция генерации значений сигнала по двум осям
        {
            double[] y = new double[length];
            PointD[] points = new PointD[length];

            for (int i = 0; i < points.Length; i++)
            {
                points[i].X = Convert.ToDouble(i);

                if (b)
                    points[i].Y = A * Math.Sin(f * 2 * Math.PI * i + Phi);
                else
                    points[i].Y = 0;
            }

            return points;
        }

        double[] GenerateNoise(PointD[] signal, double energy, bool b)  // функция генерации шума
        {
            Random rnd = new Random();
            int length = signal.Length; 

            // Генерация последовательности нормально распределённых случайных чисел.
            double[] massRand = new double[length];  // Массив случайных чисел. 
            double Er = 0;

            for (int i = 0; i < length; i++)
            {
                massRand[i] = 0;
                for (int n = 0; n < 12; n++)
                {
                    massRand[i] += Convert.ToDouble(rnd.Next(-100, 101)) / 100;  // реализация равномерно распределенной случайной величины, в интервале значений [−1, 1]
                }
                massRand[i] /= 12;  
                Er += massRand[i] * massRand[i]; 
            }

            // Подсчёт энергии шума относительно энергии сигнала.
            double Es = 0;  // Полная энергия сигнала.        
            for (int i = 0; i < signal.Length; i++)
            {
                Es += Math.Pow(signal[i].Y, 2);
            }
            double Enoise = Es * energy / 100;  // Энергия шума относительно энергии сигнала.

            // Нормировка случайной последовательности.
            for (int i = 0; i < length; i++)
            {
                if (b)
                    massRand[i] = massRand[i] * Math.Sqrt(Enoise / Er);
                else
                    massRand[i] = 0;
            }

            return massRand;
           
        }

        PointD[] SignalSummation(PointD[] signal1, PointD[] signal2, PointD[] signal3)  // функция суммирования трех сигналов
        {
            PointD[] sum = new PointD[signal1.Length];

            for (int i = 0; i < sum.Length; i++)
            {
                sum[i].X = signal1[i].X + signal2[i].X + signal3[i].X;
                sum[i].Y = signal1[i].Y + signal2[i].Y + signal3[i].Y;
            }

            return sum;
        }

        Complex[] fft(Complex[] frame, bool direct)  // БПФ
        {
            const double DoublePi = 2 * Math.PI;
            if (frame.Length == 1) return frame;
            int halfSampleSize = frame.Length / 2;
            int fullSampleSize = frame.Length;

            double arg = direct ? -DoublePi / fullSampleSize : DoublePi / fullSampleSize;
            Complex omegaPowBase = new Complex(Math.Cos(arg), Math.Sin(arg));
            Complex omega = Complex.One;
            Complex[] spectrum = new Complex[fullSampleSize];

            for (int j = 0; j < halfSampleSize; j++)
            {
                spectrum[j] = frame[j] + frame[j + halfSampleSize];
                spectrum[j + halfSampleSize] = omega * (frame[j] - frame[j + halfSampleSize]);
                omega *= omegaPowBase;
            }

            Complex[] yTop = new Complex[halfSampleSize];
            Complex[] yBottom = new Complex[halfSampleSize];
            for (int i = 0; i < halfSampleSize; i++)
            {
                yTop[i] = spectrum[i];
                yBottom[i] = spectrum[i + halfSampleSize];
            }

            yTop = fft(yTop, direct);
            yBottom = fft(yBottom, direct);
            for (int i = 0; i < halfSampleSize; i++)
            {
                int j = i << 1; 
                spectrum[j] = yTop[i];
                spectrum[j + 1] = yBottom[i];
            }

            return spectrum;
        }

        public static Complex[] FilteringSpectr(Complex[] noiseComplexSpectr, double filtrThreshold)  // очистка спектра
        {
            int lenght = noiseComplexSpectr.Length;
            Complex[] restComplexSpectr = new Complex[lenght];

            double energyThreshold = 0;
            for (int i = 0; i < lenght; i++)
                energyThreshold += noiseComplexSpectr[i].Magnitude;
            energyThreshold *= filtrThreshold / 100;

            double energy = 0;
            for (int i = 0; i < lenght >> 1; i++)
            {
                if (energy < energyThreshold)
                {
                    energy += (noiseComplexSpectr[i].Magnitude + noiseComplexSpectr[lenght - 1 - i].Magnitude);
                    restComplexSpectr[i] = noiseComplexSpectr[i];
                    restComplexSpectr[lenght - 1 - i] = noiseComplexSpectr[lenght - 1 - i];
                }
                else
                {
                    restComplexSpectr[i] = Complex.Zero;
                    restComplexSpectr[lenght - 1 - i] = Complex.Zero;
                }
            }

            return restComplexSpectr;
        }

        double StandartDeviation(PointD[] origSgnl, double[] signal_repair)  // среднеквадратичное отклонение
        {
            int lenght = origSgnl.Length;
            double SD = 0;
            for (int i = 0; i < lenght; i++)
            {
                SD = Math.Abs(signal_repair[i] - origSgnl[i].Y) * Math.Abs(signal_repair[i] - origSgnl[i].Y);
            }
            return SD / lenght;
        }

        double[] ConvertToDouble(Complex[] comp)  // преобразовать в double
        {
            double[] doub = new double[comp.Length];

            for (int i = 0; i < comp.Length; i++)
            {
                doub[i] = comp[i].Real;
            }

            return doub;
        }

        double[] ComplexMagnitude(Complex[] arr)  // функция возвращает модуль комплексных значений
        {
            double[] result = new double[arr.Length];

            for (int i = 0; i < arr.Length; i++)
                result[i] = arr[i].Magnitude;

            return result;
        }

        private void button1_Click(object sender, EventArgs e)
        {
            // Первая вкладка
            PointD[] signal_sum = SignalSummation(GeneratingPoints(Convert.ToDouble(textBox1A.Text), Convert.ToDouble(textBox1F.Text), Convert.ToDouble(textBox1P.Text), Convert.ToInt32(comboBoxFFTLength.Text), checkBox1.Checked),
                GeneratingPoints(Convert.ToDouble(textBox2A.Text), Convert.ToDouble(textBox2F.Text), Convert.ToDouble(textBox2P.Text), Convert.ToInt32(comboBoxFFTLength.Text), checkBox2.Checked),
                GeneratingPoints(Convert.ToDouble(textBox3A.Text), Convert.ToDouble(textBox3F.Text), Convert.ToDouble(textBox3P.Text), Convert.ToInt32(comboBoxFFTLength.Text), checkBox3.Checked));

            double[] noise = GenerateNoise(signal_sum, Convert.ToDouble(textBoxNoiseEnergy.Text), checkBoxNoise.Checked);

            PointD[] signal_noise_sum = new PointD[signal_sum.Length];
            for (int i = 0; i < signal_sum.Length; i++)
            {
                signal_noise_sum[i].X = Convert.ToDouble(i);
                signal_noise_sum[i].Y = noise[i] + signal_sum[i].Y;
            }

            DrawSignal(signal_noise_sum);

            // Вторая вкладка
            Complex[] signal_noise_sum_c = new Complex[signal_noise_sum.Length];
            for (int i = 0; i < signal_noise_sum_c.Length; i++)
                signal_noise_sum_c[i] = signal_noise_sum[i].Y;
            Complex[] spectrum_c = fft(signal_noise_sum_c, true);

            double[] spectrum_filtered = ComplexMagnitude(FilteringSpectr(spectrum_c, Convert.ToDouble(textBoxRepairEnergy.Text)));
            Complex[] spectrum_filtered_c = FilteringSpectr(spectrum_c, Convert.ToDouble(textBoxRepairEnergy.Text));
            DrawSpectrum(ComplexMagnitude(spectrum_c));
            DrawSpectrumClear(spectrum_filtered);

            // Третья вкладка
            Complex[] signal_repair_c = fft(spectrum_filtered_c, false);
            double[] signal_repair = ConvertToDouble(signal_repair_c);
            for (int i = 0; i < signal_repair.Length; i++)
            {
                signal_repair[i] /= signal_repair.Length;
            }
            DrawSignalResult(signal_repair);

            DrawSignalOriginal(signal_sum);

            // Вывод СКО
            textBoxCKO.Text = Convert.ToString(StandartDeviation(signal_sum, signal_repair));
        }

        private void button2_Click(object sender, EventArgs e)
        {
            Close();
        }
    }
}
