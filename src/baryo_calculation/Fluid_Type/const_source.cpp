#include <BSMPT/baryo_calculation/Fluid_Type/const_source.h>
#include <BSMPT/utility.h>

namespace BSMPT
{
    namespace Baryo
    {

        void const_source::operator()(const state_type &omega, state_type &domega, const double z)
        {
            /*
            omega[0] -> q 
            omega[1] -> t
            omega[2] -> b
            omega[3] -> l
            omega[4] -> nu
            omega[5] -> tau 
            omega[6] -> h1 
            omega[7] -> h2
            omega[8] -> u 
            omega[9]  -> q_prime
            omega[10]  -> t_prime
            omega[11]  -> b_prime
            omega[12] -> l_prime
            omega[13] -> nu_prime
            omega[14] -> tau_prime
            omega[15]  -> h1_prime
            omega[16] -> h2_prime
            omega[17] -> u_prime
        */
            /*
            Definition of all transport coefficients
       */
            std::vector<double> quark_mass;
            std::vector<double> quark_mass_prime;
            top_func(z, quark_mass, quark_mass_prime);
            double mt = quark_mass[0];
            double mb = quark_mass[1];
            std::vector<double> tau_mass, tau_mass_prime;
            tau_func(z, tau_mass, tau_mass_prime);
            double mtau = tau_mass[0];

            //Phase Calculation
            //top
            auto theta_vec_top = Calc_theta(z, gen_fluid::TOP_symmetric_CP_violating_phase, gen_fluid::TOP_broken_CP_violating_phase);
            double theta_prime_top = theta_vec_top[1];
            //bot
            auto theta_vec_bot = Calc_theta(z, gen_fluid::BOT_symmetric_CP_violating_phase, gen_fluid::BOT_broken_CP_violating_phase);
            double theta_prime_bot = theta_vec_bot[1];
            //tau
            auto theta_vec_tau = Calc_theta(z, gen_fluid::TAU_symmetric_CP_violating_phase, gen_fluid::TAU_broken_CP_violating_phase);
            double theta_prime_tau = theta_vec_tau[1];

            //TOP statistical factor
            Calc_kappa_obj.set_class(Temp, mt);
            double num_int = NIntegrate_kappa(Calc_kappa_obj);
            double kappa_tL = kappa_QL_0 * num_int;
            double kappa_tR = kappa_QR_0 * num_int;
            //BOT statistical factor
            Calc_kappa_obj.set_class(Temp, mb);
            num_int = NIntegrate_kappa(Calc_kappa_obj);
            double kappa_bL = kappa_QL_0 * num_int;
            double kappa_bR = kappa_QR_0 * num_int;
            //TAU statistical factor
            Calc_kappa_obj.set_class(Temp, mtau);
            num_int = NIntegrate_kappa(Calc_kappa_obj);
            double kappa_tauL = kappa_QL_0 * num_int;
            double kappa_tauR = kappa_QR_0 * num_int;
            //Effective statistical factor
            double kappa_q = kappa_tL * kappa_bL / (kappa_tL + kappa_bL);

            //Rescaled chemical potential like in 1811.11104
            //Relaxation Rates
            double mu_M_t = omega[1] / kappa_tR - omega[0] / kappa_q;
            double mu_M_b = omega[2] / kappa_bR - omega[0] / kappa_q;
            double mu_M_tau = omega[5] / kappa_tauR - omega[3] / kappa_tauL;
            //Yukawa Rates
            double mu_Y_t = omega[1] / kappa_tR - omega[0] / kappa_q - omega[6] / kappa_H_0 - omega[7] / kappa_H_0;
            double mu_Y_b = omega[2] / kappa_bR - omega[0] / kappa_q - omega[6] / kappa_H_0 - omega[7] / kappa_H_0;
            double mu_Y_tau = omega[5] / kappa_tauR - omega[3] / kappa_tauL + omega[6] / kappa_H_0 + omega[7] / kappa_H_0;
            //Strong sphaleron rate
            double mu_SS = -4 * omega[8] * (2 / kappa_QL_0 + 1 / kappa_QR_0) + 2 * omega[0] / kappa_q - omega[1] / kappa_tR - omega[2] / kappa_bR;
            //Numerical Integration Set up for the relaxation rates

            /*
                Additional quark-tau source term as in 1811.11104
                Eq. 15 + Eq. 14 + text below Eq.24
            */
            double mu_QL = omega[3] / kappa_tauL - omega[5] / kappa_tauR - omega[0] / kappa_q + omega[1] / kappa_tR;
            double GamQL = kappaQL * std::pow(Temp, 5) / std::pow(LambdaQL, 4);

            Calc_Gam_obj.set_class(Temp, vw, mt, msqrt_thermal_top, dmsqrt_thermal_top);
            double Gam_M_t = Nintegrate_GamM(Calc_Gam_obj);
            Calc_Gam_obj.set_class(Temp, vw, mb, msqrt_thermal_bot, dmsqrt_thermal_bot);
            double Gam_M_b = Nintegrate_GamM(Calc_Gam_obj);
            Calc_Gam_obj.set_class(Temp, vw, mtau, msqrt_thermal_tau, dmsqrt_thermal_tau);
            double Gam_M_tau = Nintegrate_GamM(Calc_Gam_obj);
            Calc_Scp_obj.set_class(Temp, vw, mt, theta_prime_top, msqrt_thermal_top, dmsqrt_thermal_top);
            double Scp_t = Nintegrate_Scp(Calc_Scp_obj);
            Calc_Scp_obj.set_class(Temp, vw, mb, theta_prime_bot, msqrt_thermal_bot, dmsqrt_thermal_bot);
            double Scp_b = Nintegrate_Scp(Calc_Scp_obj);
            Calc_Scp_obj.set_class(Temp, vw, mtau, theta_prime_tau, msqrt_thermal_tau, dmsqrt_thermal_tau);
            double Scp_tau = Nintegrate_Scp(Calc_Scp_obj);

            /*
            omega[0]    -> q 
            omega[1]    -> t
            omega[2]    -> b
            omega[3]    -> l
            omega[4]    -> nu
            omega[5]    -> tau 
            omega[6]    -> h1 
            omega[7]    -> h2
            omega[8]    -> u 
            omega[9]    -> q_prime
            omega[10]   -> t_prime
            omega[11]   -> b_prime
            omega[12]   -> l_prime
            omega[13]   -> nu_prime
            omega[14]   -> tau_prime
            omega[15]   -> h1_prime
            omega[16]   -> h2_prime
            omega[17]   -> u_prime
        */

            domega[0] = omega[9];
            domega[1] = omega[10];
            domega[2] = omega[11];
            domega[3] = omega[12];
            domega[4] = omega[13];
            domega[5] = omega[14];
            domega[6] = omega[15];
            domega[7] = omega[16];
            domega[8] = omega[17];
            /*
            dmu qmu = vw qprime - Dq q_drpime = C
            --> q_dprime = ( vw qprime - C)/Dq 
        */
            domega[9] = (vw * omega[9] - (Gam_M_t * mu_M_t + Gam_M_b * mu_M_b + Gam_Y_t * mu_Y_t + Gam_Y_b * mu_Y_b - 2 * Gam_SS * mu_SS - Scp_t - Scp_b + GamQL * mu_QL)) / Dq;
            domega[10] = (vw * omega[10] - (-Gam_M_t * mu_M_t - Gam_Y_t * mu_Y_t + Gam_SS * mu_SS + Scp_t - GamQL * mu_QL)) / Dq;
            domega[11] = (vw * omega[11] - (-Gam_M_b * mu_M_b - Gam_Y_b * mu_Y_b + Gam_SS * mu_SS + Scp_b)) / Dq;
            domega[12] = (vw * omega[12] - (Gam_M_tau * mu_M_tau + Gam_Y_tau * mu_Y_tau - Scp_tau - GamQL * mu_QL)) / Dlep;
            domega[13] = (vw * omega[13] - (0.)) / Dtau;
            domega[14] = (vw * omega[14] - (-Gam_M_tau * mu_M_tau - Gam_Y_tau * mu_Y_tau + Scp_tau + GamQL * mu_QL)) / Dtau;
            domega[15] = (vw * omega[15] - (Gam_Y_t * mu_Y_t - Gam_Y_b * mu_Y_b - Gam_Y_tau * mu_Y_tau)) / Dh;
            domega[16] = (vw * omega[16] - (Gam_Y_t * mu_Y_t - Gam_Y_b * mu_Y_b - Gam_Y_tau * mu_Y_tau)) / Dh;
            domega[17] = (vw * omega[17] - (Gam_SS * mu_SS)) / Dq;
        }
        double const_source::Calc_nL(double z_start, double z_end) const
        {
            bool debug = false;
            /*
            omega[0]    -> q 
            omega[1]    -> t
            omega[2]    -> b
            omega[3]    -> l
            omega[4]    -> nu
            omega[5]    -> tau 
            omega[6]    -> h1 
            omega[7]    -> h2
            omega[8]    -> u 
            omega[9]    -> q_prime
            omega[10]   -> t_prime
            omega[11]   -> b_prime
            omega[12]   -> l_prime
            omega[13]   -> nu_prime
            omega[14]   -> tau_prime
            omega[15]   -> h1_prime
            omega[16]   -> h2_prime
            omega[17]   -> u_prime
        */
            state_type mu(18);
            mu = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            const double C_AbsErr = 1e-9;
            const double C_RelErr = 1e-5;
            double stepsize_initial;
            if (z_start < z_end)
                stepsize_initial = 1e-8;
            if (z_start > z_end)
                stepsize_initial = -1e-8;
            double abs_err = C_AbsErr;
            double rel_err = C_RelErr;

            if (debug)
            {
                for (auto x : mu)
                    std::cout << sep << x;
                std::cout << std::endl;
            }
            integrate_adaptive(make_controlled(abs_err, rel_err, error_stepper_type()), *this, mu, z_start, z_end, stepsize_initial);
            if (debug)
            {
                for (auto x : mu)
                    std::cout << sep << x;
            }
            /*
        We have to take the sum of all left-handed quarks and leptons 
            --> q && q1 = - 2 u &&  l 
    */
            return mu[0] - 2 * mu[8] + mu[3];
        }

        void const_source::set_class(
            const int bottom_mass_inp,
            struct GSL_integration_mubl &container,
            const Calc_Gam_M &Calc_Gam_inp,
            const Calc_Scp &Calc_Scp_inp,
            const Calc_kappa_t &Calc_kappa_inp,
            const double kappaQL_inp,
            const double LambdaQL_inp)
        {

            set_class(bottom_mass_inp, container, Calc_Gam_inp, Calc_Scp_inp, Calc_kappa_inp);
            bool debug = true;
            kappaQL = kappaQL_inp;
            LambdaQL = LambdaQL_inp;
            if (debug)
            {
                std::cout << "kappaQL = " << sep << kappaQL << std::endl;
                std::cout << "LambdaQL = " << sep << LambdaQL << std::endl;
                double mt = 170;
                double theta_prime_top = TOP_broken_CP_violating_phase;
                Calc_Scp_obj.set_class(Temp, vw, mt, theta_prime_top, msqrt_thermal_top, dmsqrt_thermal_top);
                double Scp_t = Nintegrate_Scp(Calc_Scp_obj);
                std::cout << "Scp_t = " << sep << Scp_t << std::endl;
                double GamQL = kappaQL * std::pow(Temp, 5) / std::pow(LambdaQL, 4);
                std::cout << "GamQL = " << sep << GamQL << std::endl;
                std::cout << "GamYuK_t" << sep << Gam_Y_t << std::endl;
                std::cout << "GamYuk_tau" << sep <<Gam_Y_tau<<std::endl;
                Calc_Gam_obj.set_class(Temp, vw, mt, msqrt_thermal_top, dmsqrt_thermal_top);
                double Gam_M_t = Nintegrate_GamM(Calc_Gam_obj);
                std::cout << "GamM_t = "<<sep << Gam_M_t<<std::endl;
            }
        }

    } // namespace Baryo
} // namespace BSMPT