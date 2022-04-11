using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;


namespace MyApp
{
    class Program
    {
        static void Main(string[] args)
        {
            Exchange e = new Exchange();
            Console.Write("Please Enter the Exchange Name: ");
            e.Name = Console.ReadLine();
            Console.Write("Please Enter the Exchange Ticker: ");
            e.Symbol = Console.ReadLine();

            Underlying u = new Underlying();
            u.Exchange = e;
            Console.Write("What is the Underlying Ticker: ");
            u.Symbol = Console.ReadLine();
            Console.Write("What is the Underlying Price: ");
            u.Price = Convert.ToDouble(Console.ReadLine());


            EuropeanOption euro = new EuropeanOption();
            Console.Write("Is Call? (true or false): ");
            euro.IsCall = Convert.ToBoolean(Console.ReadLine());
            Console.Write("What is the Strike Price: ");
            euro.Strike = Convert.ToDouble(Console.ReadLine());

            Console.Write("What year the Option Expire: ");
            int Year = Convert.ToInt32(Console.ReadLine());

            Console.Write("What Month the Option Expire: ");
            int Month = Convert.ToInt32(Console.ReadLine());

            Console.Write("What Day the Option Expire: ");
            int Day = Convert.ToInt32(Console.ReadLine());

            euro.ExpirationDate = new DateTime(Year, Month, Day);
            euro.Underlying = u;

            Volatility vol = new Volatility();
            Console.Write("What is the volatility: ");
            vol.Vol = Convert.ToDouble(Console.ReadLine());

            Console.Write("What is the risk-free interest rate: ");
            double rate = Convert.ToDouble(Console.ReadLine());

            Console.Write("How many steps you would like to simulate: ");
            int Steps = Convert.ToInt32(Console.ReadLine());

            Console.Write("How many sitmulations you would like to simulate: ");
            int Simulations = Convert.ToInt32(Console.ReadLine());



            var g = euro.GetPriceAndGreeks(Steps, Simulations, vol, rate);
            if (euro.IsCall == true)
            {
                Console.WriteLine("For " + u.Symbol +" European Call option which expire on " +euro.ExpirationDate+":");
            }
            else if (euro.IsCall == false)
            {
                Console.WriteLine("For " + u.Symbol +" European Put option which expire on " +euro.ExpirationDate+":");
            }
            Console.WriteLine("Option Price: ");
            Console.WriteLine(g.Price);
            Console.WriteLine("Option Delta: ");
            Console.WriteLine(g.Delta);
            Console.WriteLine("Option Gamma: ");
            Console.WriteLine(g.Gamma);
            Console.WriteLine("Option Vega: ");
            Console.WriteLine(g.Vega);
            Console.WriteLine("Option Theta: ");
            Console.WriteLine(g.Theta);
            Console.WriteLine("Option Rho: ");
            Console.WriteLine(g.Rho);
        }
    }
    class RatePoints
    {
        public double Tenor {get; set;}
        public double Rate {get; set;}
    }
    class YeiildCurve
    {
        public List<RatePoints> CurvePoints {get; set;}
    }
    class Exchange
    {
        public string Name{ get; set;}
        public string Symbol{get; set;}
    }
    class Underlying
    {
        public string Symbol {get; set;}
        public Exchange Exchange{get; set;}
        public double Price {get; set;}
    }
    abstract class Option
    {
        public DateTime ExpirationDate {get; set;}
        public Underlying Underlying {get; set;}
        public abstract OptionResult GetPriceAndGreeks  (long Steps, long Simulations, Volatility vol, double rate);
    }
    class OptionResult 
    {
        public double Price {get; set;}
        public double Delta{get; set;}
        public double Theta{get; set;}
        public double Gamma{get; set;}
        public double Vega{get; set;}
        public double Rho{get; set;}
    }
    class Volatility
    {
        public double Vol{get; set;}
    }
    
    class EuropeanOption : Option
    {
        public double Strike {get; set;}
        public bool IsCall {get; set;}

        public override OptionResult GetPriceAndGreeks(long Steps, long Simulations, Volatility vol, double rate)
        {
            GaussianRandoms r = new GaussianRandoms();
            r.PopulatedNRand(122, Simulations, Steps);

            SimulationParameters p = new SimulationParameters();
            p.S0 = Underlying.Price;
            p.r = rate; // Yeild curve function
            p.Steps = Steps;
            p.Simulations = Simulations;
            p.Tenor = (ExpirationDate - DateTime.Today).Days / 365d;
            
            p.Volatility = vol.Vol;
            OptionResult PriceAndGreeks = new OptionResult();
            
            
            static double OptionPrice(SimulationParameters p, GaussianRandoms r,double Strike, bool IsCall, long Simulations, long Steps) 
            {   

                var result = MonteCarloSimulator.GeneratePaths(p,r);
                double discountFactor = Math.Exp(-p.r*p.Tenor);
                double STsum = 0;
                if (IsCall == true)
                {
                    for (int i = 0; i < Simulations; i++)
                    {
                        STsum = STsum + Math.Max(result.SimulatedPaths[i, Steps - 1]- Strike,0);
                    
                    }
                }
                else if (IsCall == false)
                {
                    for (int i = 0; i < Simulations; i++)
                    {
                        STsum = STsum + Math.Max(Strike - result.SimulatedPaths[i, Steps - 1],0);
                    
                    }
                }
                return (discountFactor*STsum/Simulations);
            }
            PriceAndGreeks.Price = OptionPrice(p,r,Strike,IsCall, Simulations, Steps);

            double deltaS = 0.001;
            double deltaSigma = 0.01;
            double deltaT = 0.00273972602;
            double deltar = 0.01;
            SimulationParameters p_deltaSPlus = new SimulationParameters();
            p_deltaSPlus = p;
            p_deltaSPlus.S0 = p.S0 + deltaS;
            double OptionPriceDeltaPlus = OptionPrice(p_deltaSPlus,r,Strike,IsCall,Simulations,Steps);
            SimulationParameters p_deltaSMinus = new SimulationParameters();
            p_deltaSMinus = p;
            p_deltaSMinus.S0 = p.S0 - deltaS;
            double OptionPriceDeltaMinus = OptionPrice(p_deltaSMinus,r,Strike,IsCall,Simulations,Steps);
            PriceAndGreeks.Delta = (OptionPriceDeltaPlus - OptionPriceDeltaMinus)/(2*deltaS);
            PriceAndGreeks.Gamma = (OptionPriceDeltaPlus - 2*PriceAndGreeks.Price + OptionPriceDeltaMinus)/(Math.Pow(deltaS,2));

            SimulationParameters p_deltaVolPlus = new SimulationParameters();
            p_deltaVolPlus = p;
            p_deltaVolPlus.Volatility = p.Volatility + deltaSigma;
            double OptionPriceVolPlus = OptionPrice(p_deltaVolPlus,r,Strike,IsCall,Simulations,Steps);
            SimulationParameters p_deltaVolMinus = new SimulationParameters();
            p_deltaVolMinus = p;
            p_deltaVolMinus.Volatility = p.Volatility - deltaSigma;
            double OptionPriceVolMinus = OptionPrice(p_deltaVolMinus,r,Strike,IsCall,Simulations,Steps);
            PriceAndGreeks.Vega = (OptionPriceVolPlus - OptionPriceVolMinus)/(2 * deltaSigma);

            SimulationParameters p_deltatPlus = new SimulationParameters();
            p_deltatPlus = p;
            p_deltatPlus.Volatility = p.Tenor + deltaT;
            double OptionPricetPlus = OptionPrice(p_deltatPlus,r,Strike,IsCall,Simulations,Steps);
            PriceAndGreeks.Theta = (OptionPricetPlus - PriceAndGreeks.Price)/deltaT;

            SimulationParameters p_deltarPlus = new SimulationParameters();
            p_deltarPlus = p;
            p_deltarPlus.r = p.r + deltar;
            double OptionPricerPlus = OptionPrice(p_deltarPlus,r,Strike,IsCall,Simulations,Steps);
            SimulationParameters p_deltarMinus = new SimulationParameters();
            p_deltarMinus = p;
            p_deltarMinus.Volatility = p.r - deltar;
            double OptionPricerMinus = OptionPrice(p_deltarMinus,r,Strike,IsCall,Simulations,Steps);
            PriceAndGreeks.Rho = (OptionPricerPlus - OptionPricerMinus)/(2 * deltar);


            return PriceAndGreeks;
        }
    }
    

    class GaussianRandoms
    {
        public double[,] NRands {get; set;}
        public void PopulatedNRand(int seed, long rows, long cols)
        {
            NRands = new double [rows, cols];
            Random r_1 = new Random(seed);
            Random r_2 = new Random(seed);
            for(int i = 0; i <rows; i++)
            {
                for (int j = 0; j< cols; j++)
                {
                    double x1 = r_1.NextDouble();
                    double x2 = r_2.NextDouble();
                    double z1 = Math.Sqrt(-2 * Math.Log(x1))*Math.Cos(2*Math.PI*x2);
                    double z2 = Math.Sqrt(-2 * Math.Log(x1))*Math.Sin(2*Math.PI*x2);
                    NRands[i,j] = z1;
                }
            }
        }
    }
    class SimulationResult
    {
        public double[,] SimulatedPaths {get; set;}
    }
    class SimulationParameters
    {
        public double S0{get; set;}
        public double r{get; set;}
        public long Steps{get; set;}
        public long Simulations{get; set;}
        public double Tenor {get; set;}
        public double Volatility{get; set;}
    }
    static class MonteCarloSimulator
    {
        public static SimulationResult GeneratePaths(SimulationParameters p, GaussianRandoms rands)
        {
            
            SimulationResult results = new SimulationResult();
            results.SimulatedPaths = new double [p.Simulations,p.Steps];
            double dt = p.Tenor/p.Steps;
            
            for (int i = 0; i < p.Simulations; i++)
            {
                for (int j = 0; j < p.Steps; j++)
                {
                    if (j == 0)
                    {
                        results.SimulatedPaths[i,j] = p.S0;
                    }
                    else
                    {
                        results.SimulatedPaths[i,j] = results.SimulatedPaths[i, j-1] * Math.Exp(((p.r - 0.5 * Math.Pow(p.Volatility,2))* dt) + (p.Volatility * Math.Sqrt(dt)*rands.NRands[i,j]));
                    }
                }
            }
            return results;
        }
    }
}
//Colab with Hu Hong
