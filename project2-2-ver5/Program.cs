using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Numerics;
using static System.Console;
using static System.Math;

namespace project2_2_ver5
{
    class Program
    {
        //
        // 두께에 대한 MSE를 구하기 위한 클래스
        //
        // 2021.04.01 정지훈.
        //
        class MSE
        {
            private int LenData;

            private double[] n_SiO2;
            private double[] k_SiO2;

            // 반사계수를 담을 배열.
            private Complex[] r12p;
            private Complex[] r12s;
            private Complex[] r01p;
            private Complex[] r01s;

            // 투과계수를 담을 배열.
            private Complex[] t12p;
            private Complex[] t12s;
            private Complex[] t01p;
            private Complex[] t01s;

            private double AOI_air;
            private Complex N_air;

            private double[] wavelength_SiO2;

            private List<double> alpha_exp = new List<double>();
            private List<double> beta_exp = new List<double>();

            public MSE()
            {

                List<string> MeasurementSpectrumData = new List<string>();  // 측정 스펙트럼 데이터 저장할 배열. (한 줄씩 저장)
                string[] SingleLineData;                                    // 한 줄의 스펙트럼 데이터를 임시로 저장할 배열.

                // "SiO2 2nm_on_Si.dat" 파일 읽기. (한 줄씩)
                MeasurementSpectrumData.AddRange(File.ReadAllLines("SiO2 1000nm_on_Si.dat"));

                // 무의미한 공백 행을 제거한다.
                int lenSpectrumData = MeasurementSpectrumData.Count;
                string Blank = "";
                for (int i = 0; i < lenSpectrumData; i++)
                {
                    if (MeasurementSpectrumData[0] == Blank)
                        MeasurementSpectrumData.RemoveAt(0);
                    else
                        break;
                }

                // wavelength : 350 ~ 980(nm)인 측정 스펙트럼 데이터를 담을 리스트 선언.
                List<double> wavelength_exp = new List<double>();   // 파장 데이터 리스트.
                List<double> AOI_exp = new List<double>();          // 입사각 데이터 리스트.
                List<double> alpha_exp = new List<double>();        // Psi 데이터 리스트.
                List<double> beta_exp = new List<double>();         // Delta 데이터 리스트.

                // 데이터의 첫번째 줄은 column 명이다.
                // 이를 제외하기 위해 반복문을 1부터 시작한다.
                int StartIndex = 1;
                int LenData = MeasurementSpectrumData.Count;
                for (int i = StartIndex; i < LenData; i++)
                {
                    // tsv 형식의 데이터를 SingleLineData에 저장한다.
                    SingleLineData = MeasurementSpectrumData[i].Split((char)0x09);  // 0x09 : 수평 탭.

                    // 파장이 350 ~ 980(nm) 이내인 데이터만 저장한다.
                    if (Convert.ToDouble(SingleLineData[0]) >= 350.0 &&
                        Convert.ToDouble(SingleLineData[0]) <= 980.0)
                    {
                        // 각 컬럼에 해당하는 데이터를 저장한다.
                        wavelength_exp.Add(Double.Parse(SingleLineData[0]));
                        AOI_exp.Add(Double.Parse(SingleLineData[1]));
                        alpha_exp.Add(Double.Parse(SingleLineData[2]));
                        beta_exp.Add(Double.Parse(SingleLineData[3]));
                    }
                }

                // psi, delta -> alpha, beta 변환.
                // degree, radian 변환 인라인 함수 정의.
                double degree2radian(double angle) => ((angle * (PI)) / 180.0);

                // Polarizer offset 각도. (45도)
                double PolarizerRadian = degree2radian(45.0);

                // psi, delta 데이터를 alpha, beta 로 변환한다.
                LenData = wavelength_exp.Count;


                for (int i = 0; i < LenData; i++)
                {
                    // psi, delta 값을 radian 으로 변환한다.
                    double PsiRadian = degree2radian(alpha_exp[i]);
                    double DeltaRadian = degree2radian(beta_exp[i]);

                    // psi, delta 데이터를 alpha, beta 로 갱신한다.
                    alpha_exp[i] = (
                        (Pow(Tan(PsiRadian), 2.0) - Pow(Tan(PolarizerRadian), 2.0))
                        / (Pow(Tan(PsiRadian), 2.0) + Pow(Tan(PolarizerRadian), 2.0)));
                    beta_exp[i] = (
                        (2.0 * Tan(PsiRadian) * Tan(PolarizerRadian) * Cos(DeltaRadian))
                        / (Pow(Tan(PsiRadian), 2.0) + Pow(Tan(PolarizerRadian), 2.0)));
                    this.alpha_exp.Add(alpha_exp[i]);
                    this.beta_exp.Add(beta_exp[i]);
                }

                // "Si_new.txt" 파일 읽기.
                string[] Si_new = File.ReadAllLines("Si_new.txt");  // Si 기판 물성값 저장.(한 줄씩)

                // 데이터의 첫번째 줄은 column 명이다.
                // 이를 제외하고 데이터를 받기 위해 LenData 변수를 선언한다.
                LenData = Si_new.Length - 1;
                double[] wavelength_Si = new double[LenData];
                double[] n_Si = new double[LenData];
                double[] k_Si = new double[LenData];

                // Si_new 에 받은 데이터를 각 컬럼별로 저장한다.
                LenData = Si_new.Length;

                for (int i = StartIndex; i < LenData; i++)
                {
                    // tsv 형식의 데이터를 SingleLineData에 저장한다.
                    SingleLineData = Si_new[i].Split((char)0x09);  // 0x09 : 수평 탭.

                    // 각 컬럼에 해당하는 데이터를 저장한다.
                    wavelength_Si[i - 1] = Double.Parse(SingleLineData[0]);
                    n_Si[i - 1] = Double.Parse(SingleLineData[1]);
                    k_Si[i - 1] = Double.Parse(SingleLineData[2]);
                }


                // "SiO2_new.txt" 파일 읽기.
                string[] SiO2_new = File.ReadAllLines("SiO2_new.txt");  // Si 기판 물성값 저장.(한 줄씩)

                // 데이터의 첫번째 줄은 column 명이다.
                // 이를 제외하고 데이터를 받기 위해 LenData 변수를 선언한다.
                LenData = SiO2_new.Length - 1;
                double[] wavelength_SiO2 = new double[LenData];
                double[] n_SiO2 = new double[LenData];
                double[] k_SiO2 = new double[LenData];

                this.n_SiO2 = new double[LenData];
                this.k_SiO2 = new double[LenData];
                this.wavelength_SiO2 = new double[LenData];

                // SiO2_new 에 받은 데이터를 각 컬럼별로 저장한다.
                LenData = SiO2_new.Length;
                for (int i = StartIndex; i < LenData; i++)
                {
                    // tsv 형식의 데이터를 SingleLineData에 저장한다.
                    SingleLineData = SiO2_new[i].Split((char)0x09);  // 0x09 : 수평 탭.

                    // 각 컬럼에 해당하는 데이터를 저장한다.
                    wavelength_SiO2[i - 1] = Double.Parse(SingleLineData[0]);
                    n_SiO2[i - 1] = Double.Parse(SingleLineData[1]);
                    k_SiO2[i - 1] = Double.Parse(SingleLineData[2]);

                    this.n_SiO2[i - 1] = n_SiO2[i - 1];
                    this.k_SiO2[i - 1] = k_SiO2[i - 1];
                    this.wavelength_SiO2[i - 1] = wavelength_SiO2[i - 1];
                }

                LenData = wavelength_Si.Length;

                // 반사계수를 담을 배열.
                Complex[] r12p = new Complex[LenData],
                          r12s = new Complex[LenData],
                          r01p = new Complex[LenData],
                          r01s = new Complex[LenData];
                // 투과계수를 담을 배열.
                Complex[] t12p = new Complex[LenData],
                          t12s = new Complex[LenData],
                          t01p = new Complex[LenData],
                          t01s = new Complex[LenData];

                // 반사계수를 담을 배열.
                this.r12p = new Complex[LenData];
                this.r12s = new Complex[LenData];
                this.r01p = new Complex[LenData];
                this.r01s = new Complex[LenData];
                // 투과계수를 담을 배열.
                this.t12p = new Complex[LenData];
                this.t12s = new Complex[LenData];
                this.t01p = new Complex[LenData];
                this.t01s = new Complex[LenData];

                double AOI_air = degree2radian(65.0);   // 입사각. (라디안) 
                Complex N_air = new Complex(1.0, 0);    // 공기의 굴절률.
                this.AOI_air = AOI_air;
                this.N_air = N_air;
                // 반사, 투과계수를 계산한다.
                for (int i = 0; i < LenData; i++)
                {
                    // 파장에 대한 물질의 복소굴절률을 구한다.
                    Complex N_SiO2 = new Complex(n_SiO2[i], -k_SiO2[i]);
                    Complex N_Si = new Complex(n_Si[i], -k_Si[i]);

                    // air, SiO2 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                    Complex Sintheta_j = new Complex(Sin(AOI_air), 0);
                    Complex Costheta_j = new Complex(Cos(AOI_air), 0);
                    Complex Sintheta_k = (N_air / N_SiO2) * Sintheta_j;
                    Complex theta_k = Complex.Asin(Sintheta_k);
                    // air, SiO2 경계면에서의 굴절각.
                    Complex Costheta_k = Complex.Cos(theta_k);

                    // air, SiO2 경계면에서의 반사계수를 구한다.
                    r01p[i] = ((N_SiO2 * Costheta_j) - (N_air * Costheta_k)) /
                                   ((N_SiO2 * Costheta_j) + (N_air * Costheta_k));

                    r01s[i] = ((N_air * Costheta_j) - (N_SiO2 * Costheta_k)) /
                                   ((N_air * Costheta_j) + (N_SiO2 * Costheta_k));
                    this.r01p[i] = r01p[i];
                    this.r01s[i] = r01s[i];
                    // air, SiO2 경계면에서의 투과계수를 구한다.
                    t01p[i] = (N_air * Costheta_j * 2.0) /
                                   ((N_SiO2 * Costheta_j) + (N_air * Costheta_k));

                    t01s[i] = (N_air * Costheta_j * 2.0) /
                                   ((N_air * Costheta_j) + (N_SiO2 * Costheta_k));
                    this.t01p[i] = t01p[i];
                    this.t01s[i] = t01s[i];
                    // SiO2, Si 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                    Sintheta_j = Complex.Sin(theta_k);
                    Costheta_j = Complex.Cos(theta_k);
                    Sintheta_k = (N_SiO2 / N_Si) * Sintheta_j;
                    theta_k = Complex.Asin(Sintheta_k);         // SiO2, Si 경계면에서의 굴절각.
                    Costheta_k = Complex.Cos(theta_k);

                    // SiO2, Si 경계면에서의 반사계수를 구한다.
                    r12p[i] = ((N_Si * Costheta_j) - (N_SiO2 * Costheta_k)) /
                                 ((N_Si * Costheta_j) + (N_SiO2 * Costheta_k));

                    r12s[i] = ((N_SiO2 * Costheta_j) - (N_Si * Costheta_k)) /
                                 ((N_SiO2 * Costheta_j) + (N_Si * Costheta_k));
                    this.r12p[i] = r12p[i];
                    this.r12s[i] = r12s[i];
                    // SiO2, Si 경계면에서의 투과계수를 구한다.
                    t12p[i] = (N_SiO2 * Costheta_j * 2.0) /
                                 ((N_Si * Costheta_j) + (N_SiO2 * Costheta_k));

                    t12s[i] = (N_SiO2 * Costheta_j * 2.0) /
                                 ((N_SiO2 * Costheta_j) + (N_Si * Costheta_k));
                    this.t12p[i] = t12p[i];
                    this.t12s[i] = t12s[i];
                }
                this.LenData = LenData;
            }

            public double returnMse(double thickness)
            {


                // 두께 범위와 두께 간격을 설정한다.
                double gap = 5.0;

                // MSE 와 두께를 담을 배열을 선언, 초기화한다.
                double MSEs;

                // 두께별 MSE 를 계산해서 MSEs 배열에 저장한다.

                // 총 반사계수를 저장할 배열 선언.
                Complex[] Rp = new Complex[LenData],
                          Rs = new Complex[LenData];

                for (int i = 0; i < LenData; i++)
                {
                    // SiO2의 복소 굴절률.
                    Complex N_SiO2 = new Complex(n_SiO2[i], -k_SiO2[i]);

                    // air, SiO2 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                    Complex Sintheta_j = new Complex(Sin((double)AOI_air), 0);
                    Complex Sintheta_k = (N_air / N_SiO2) * Sintheta_j;
                    Complex theta_k = Complex.Asin(Sintheta_k);         // air, SiO2 경계면에서의 굴절각.
                    Complex Costheta_k = Complex.Cos(theta_k);

                    // 위상 두께를 구한다.
                    Complex PhaseThickness = ((double)thickness * Math.PI * 2.0 / wavelength_SiO2[i]) * N_SiO2 * Costheta_k;

                    //WriteLine(PhaseThickness);

                    // 총 반사계수를 구한다.
                    Complex E = Complex.Exp(PhaseThickness * new Complex(0, -2.0));

                    Rp[i] = (r01p[i] + r12p[i] * E) /
                            (1 + r01p[i] * r12p[i] * E);

                    Rs[i] = (r01s[i] + r12s[i] * E) /
                            (1 + r01s[i] * r12s[i] * E);
                }

                // alpha, beta 이론값을 담을 배열 선언.
                double[] alpha_cal = new double[LenData],
                         beta_cal = new double[LenData];

                // Polarizer 오프셋 각.
                double degree2radian(double angle) => ((angle * (PI)) / 180.0);
                double polarizerAngle = degree2radian(45.0);

                for (int i = 0; i < LenData; i++)
                {
                    // 총 반사계수비. (복소반사계수비)
                    Complex rho = Rp[i] / Rs[i];

                    // Psi, Delta.
                    double Psi = Atan(rho.Magnitude);
                    double Delta = rho.Phase;


                    alpha_cal[i] = (Pow(Tan(Psi), 2.0) - Pow(Tan(polarizerAngle), 2.0)) /
                                           (Pow(Tan(Psi), 2.0) + Pow(Tan(polarizerAngle), 2.0));

                    beta_cal[i] = (2.0 * Tan(Psi) * Cos(Delta) * Tan(polarizerAngle)) /
                                           (Pow(Tan(Psi), 2.0) + Pow(Tan(polarizerAngle), 2.0));
                }
                double sum = 0;
                for (int i = 0; i < LenData; i++)
                {
                    double difference_MSE =
                         Pow((alpha_exp[i] - alpha_cal[i]), 2.0) +
                         Pow((beta_exp[i] - beta_cal[i]), 2.0);
                    sum += difference_MSE;

                }

                MSEs = sum / LenData;
                //WriteLine($"{ thickness }     { MSEs }");
                return MSEs;
            }
        }

        static void Main(string[] args)
        {
            MSE MyMSE = new MSE();

            double StartThickness = 700.0;
            double EndThickness = 1300.0;
            double gap = 5.0;

            // MSE 와 두께를 담을 배열을 선언, 초기화한다.
            double numMSE = (EndThickness - StartThickness) / gap + 1;
            double[] MSEs = new double[(int)numMSE];
            double[] thicknesses = new double[(int)numMSE];

            // 두께별 MSE 를 계산해서 MSEs 배열에 저장한다.
            int idx = 0;
            for (double thickness = StartThickness; thickness <= EndThickness; thickness += gap)
            {
                thicknesses[idx] = thickness;
                MSEs[idx] = MyMSE.returnMse(thickness);
                ++idx;
            }

            //
            // 두께별로 계산된 MSE 를 통해
            // global minimum 일 때의 두께, MSE 값을 찾는다.
            //
            // 2021.03.26 이지원.
            //
            #region global minimum 탐색.

            // MSEs 에서 global minimum 에서의 MSE, 두께값을 구한다.
            int idxGlobalMinimum = 0;           // global minimum 에서의 index.
            double GlobalMinimum = MSEs.Min();  // global minimum 값.
            int LenData;

            // global minimum 에서의 index 를 찾는다.
            LenData = MSEs.Length;
            for (int i = 0; i < LenData; i++)
            {
                if (MSEs[i] == GlobalMinimum)
                {
                    idxGlobalMinimum = i;
                    break;
                }
            }

            //WriteLine($"{thicknesses[idxGlobalMinimum]}     {MSEs[idxGlobalMinimum]}     {idxGlobalMinimum}");
            //1070     0.0016568615951739671     74

            // global minimum 에서의 두께 값.
            double nowD0 = thicknesses[idxGlobalMinimum];
            double globalD0 = thicknesses[idxGlobalMinimum];
            double globalMSE = MSEs[idxGlobalMinimum];

            #endregion

            //
            // d0 근처에 있는 "실제" global minimum 이 되는 두께 d_sol 을 찾는다.
            // MSE 값의 변화가 10^-5 이하가 되면, 최적화 수행을 멈춘다.
            //
            // 2021.03.29 이지원.
            //
            #region 두께 범위를 재설정하고 MSE를 구한다.
            double d_sol = 0.0;             // "실제" global minimum 이 되는 두께.

            // 이전 값을 가져야하는 변수에 현재 값을 초기화 해둔다.
            double preD0 = nowD0,           // 이전 두께를 저장하는 변수.
                   tempD0 = 0.0;
            double preMSE = GlobalMinimum,  // 이전 MSE 를 저장하는 변수.
                   nowMSE = GlobalMinimum;  // 현재 MSE 를 저장하는 변수.
            double preMSE_gradient = GlobalMinimum,  // 이전 MSE 를 저장하는 변수.
                   nowMSE_gradient = GlobalMinimum;  // 현재 MSE 를 저장하는 변수.

            // global minimum 을 찾는 시간을 측정한다.
            Stopwatch stopwatch = new Stopwatch();

            ///////////////////////////////////////////////////////////////////////////////////////////////////
            nowD0 = 800;
            preD0 = 0;

            // 두께 간격, 현재 두께를 갱신한다.
            gap = 0.5;

            double gt = 0; // 
            double mt = 0; // 아담
            double[] gradient = new double[2] { 0, 0 }; // 기울기
            int cnt = 0;    // global minimum 을 찾기 위한 연산의 반복 횟수.
            double vt = 0;

            MSE mSE = new MSE();
            preMSE = mSE.returnMse(0);

            stopwatch.Start();
            while (true)
            {
                ++cnt;
                WriteLine($"==== {cnt} 회차 ====");

                // 현재 두께의 mse를 구한다.
                nowMSE = mSE.returnMse(nowD0);

                WriteLine($"predo : {preD0}     preMSE: {preMSE}\n" +
                    $"nowdo : {nowD0}     nowMSE: {nowMSE}\n" +
                    $"beforegradient : {gradient[0]}     gradient : {gradient[1]}\n" +
                    $"gap : {gap}     vt : {vt}");

                bool a = true, b = true;
                _ = (gradient[0] >= 0) ? a = false : b = false;
                _ = (gradient[1] >= 0) ? a = false : b = false;

                // 두 MSE 의 차이가 0.00001(10^-5) 이하이면 while 문을 탈출한다.
                if (Abs(preMSE - nowMSE) <= 0.00001 && a == false && b == false && preD0 != nowD0)
                {
                    // global minimum 은 둘 중 작은 MSE 로 정한다.
                    GlobalMinimum = (preMSE < nowMSE) ? preMSE : nowMSE;
                    // 실제 두께는 global minimum 으로 정해진 MSE 와 짝이 되는 두께이다.
                    d_sol = (preMSE == GlobalMinimum) ? preD0 : nowD0;
                    goto Find_d_sol;
                }

                gradient[0] = gradient[1];
                gradient[1] = (globalMSE - nowMSE) / (globalD0 - nowD0);
                preMSE = nowMSE;
                preD0 = nowD0;

                //아담 수정2
                //
                gt = gradient[0];
                mt = 0.9 * mt + (1 - 0.9) * gt;
                vt = 0.999 * vt + (1 - 0.999) * Math.Pow(gradient[0], 2);
                double mt_cap = mt / (1 - Math.Pow(0.9, cnt));
                double vt_cap = vt / (1 - Math.Pow(0.999, cnt));
                nowD0 = preD0 - (gap * mt_cap) / (Sqrt(vt_cap) + Math.Pow(10, -8));

                //베타2 = 0.999 베타1 = 0.9 알파 초깃값 = 0.001
                //gap = (0.001 * Sqrt(1-Math.Pow(0.999, cnt))) / (1- Math.Pow(0.9, cnt));
            }

        Find_d_sol: WriteLine($"\n\n\n d_sol: {d_sol}    MSE: {GlobalMinimum}    idxGlobalMinimum: {idxGlobalMinimum}");
            stopwatch.Stop();
            WriteLine($"소요 시간: {stopwatch.ElapsedMilliseconds}ms");

            #endregion
        }
    }
}
