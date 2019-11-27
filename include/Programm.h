#include <cmath>
#include <iostream>
#include <utility>
#include <stdexcept>
#include <vector>
#include <random>
#include <algorithm>

#include "Save.h"
/*
	Структура для описания условий задачи
*/
struct Data
	{
		std::pair<double,double> X_interval = {0.0, 3.141592653589793238463}; // Интервал
		size_t K_ = 100; // Количество отсчетов.
		double amplitude_of_noise = 0.5; // Амплитуда равномерного шума
		double P = 0.95; // Вероятность попадания в окрестность экстремума
		double eps = 0.01; // Интервал неопределенности
		double J = 90000000.0;
		size_t L_ = 10u; // Интервал неопределенности
		std::vector <size_t> radius_ {3,5}; // Размер скользящего окна
		std::vector <size_t> M_;
		std::vector <double> alpha;
		std::vector <double> xk;
		std::vector <double> F_xk_;
		std::vector <double> F_xk_noised;
		std::vector <double> F_xk_filtered;
		std::vector <double> lambda; // Веса свертки
		std::pair<double,double> omega_delta{90000000.0,90000000.0};

		void init()
			{
				std::random_device rd;
				std::mt19937 gen(rd());
				std::uniform_real_distribution <double> dis (-amplitude_of_noise/2.0,amplitude_of_noise/2.0);

				for (size_t index = 0u;index <= L_;++index)
					{
						lambda.push_back(0.1 * index);
					}
				// Вычисляем xk
				for (size_t index = 0u;index <= K_;++index)
					{
						xk.push_back(calculate_x_k(index));
					}
				// Вычисляем f(xk) и f(xk) с шумами
				for (size_t index = 0u;index <= K_;++index)
					{
						F_xk_.push_back(signal(xk[index]));
						F_xk_noised.push_back(signal_plus_noise(xk[index],dis(gen)));
					}
				// Вычисляем M
				for (size_t index = 0u;index < radius_.size(); ++index)
					{
						M_.push_back((radius_[index]-1u)/2u);
					}
				std::cout<<"";
			}

		const double calculate_x_k(size_t k)
			{
				if (k >= 0u && k <= 100u)
					{
						return X_interval.first + k * (X_interval.second - X_interval.first)/K_;
					}
				else
					{
						throw std::runtime_error("Неверное значение для k!");
					}
			}
		/**
			Функция вычисления сигнала для xk
		*/
		const double signal(double Xk_)
			{
				return std::sin(Xk_) + 0.5;
			}
		/**
			Функция вычисления сигнала с шумом для xk
		*/
		const double signal_plus_noise(double Xk_,double noise)
			{
				return std::sin(Xk_) + 0.5 + noise;
			}
		/**
			Функция получения числа испытаний N 
		*/
		const double get_N()
			{
				return (log(1.0-P)/log(1.0 - ( eps/( X_interval.second - X_interval.first ) )));
			}
		/**
			Функция вычисления среднего геометрического
		*/
		double geom_average(std::vector<double> alpha,size_t M,size_t k)
			{
				double temp = 1.0;
				if(k < M || k > K_ - M)	// k-M должно быть больше или равно 0
					{
						return 0.0;
					}
				size_t degree = 0u;

				for (size_t index = k - M ; index <= k + M ;++index)
					{
						degree = index + M - k ;
						if (degree >= alpha.size() )
							{
								degree  = alpha[ alpha.size()-(degree-(alpha.size()-1u))-1u ];
							}
						temp *= pow(F_xk_noised[index],alpha[degree]);
					}
				return temp;
			}
		double sum(std::vector<double> alpha,size_t a,size_t b)
			{
				double temp = 0.0;
				for (size_t index = a;index <= b;++index)
					{
						if (index + 1u > alpha.size())
							{
								temp+=alpha[index-alpha.size()];
								continue;
							}
						temp+=alpha[index];
					}
				return temp;
			}
		void init_Fx_filtered(std::vector <double> input_alpha,size_t M)
			{
				F_xk_filtered.clear();
				for (size_t index = 0u;index <= K_;++index)
					{
						F_xk_filtered.push_back(geom_average(input_alpha,M,index));
					}
			}
	};

class worker
	{
		double N_;
		V10 v10;
		std::vector <Save> save;
		Data data;
	public:
		worker()
			{
				data.init();
				N_ = data.get_N();	
			}
		void pass()
			{
				Save temp_save;
				std::pair<double,double> temp_omega_delta;
				double temp_J = 0.0;
				std::random_device rd;
				std::mt19937 gen(rd());
				std::uniform_real_distribution <double> dis (0.0, 1.0);
				std::vector<double> temp_alpha;

				double temp;

				temp_save.graphic.x = data.xk;
				temp_save.graphic.Function = data.F_xk_;
				temp_save.graphic.Function_noised = data.F_xk_noised;
				
				for (size_t jump = 0u; jump < data.radius_.size(); ++jump)
					{
						double r = data.radius_[jump];
						double M = data.M_[jump];

						temp_save.r = data.radius_[jump];
						for (size_t l = 0u; l < data.lambda.size(); ++l) // для различных значений весов
							{
								temp_save.h = data.lambda[l];
								for (size_t index = 0u; index < N_; ++index)		// Число испытаний N
									{
										temp_alpha.reserve(M+1u);
										temp_alpha.resize(M+1u);
										temp_alpha.back() = dis(gen);
										if (data.M_[jump] >= 2u)
											{
												for (size_t count = 1u; count < data.M_[jump] ; ++count)
													{
														std::uniform_real_distribution <double> temp_dis(0.0,1.0 - data.sum(temp_alpha,M, r-M-1));
														temp_alpha.at(count) = 0.5 * temp_dis(gen);
													}
											}
										temp_alpha.front() = 0.5*(1.0 - data.sum(temp_alpha,1u,r-2u));
										// Уже сгенерирован набор альфа
										data.init_Fx_filtered(temp_alpha,M);	// Находим отфильтрованную функцию
										// Находим критерий зашумленности (по Евклиду)
										double summ = 0.0;
										for (size_t kol = 1u; kol <= data.K_; ++kol) 
											{
												summ += pow((data.F_xk_filtered[kol] - data.F_xk_filtered[kol - 1u]) , 2);		
											}
										temp_omega_delta.first = sqrt(summ); summ=0.0;
										// Находим критерий отличия (по Евклиду)
										for (size_t kol = 0u; kol <= data.K_; ++kol) 
											{
												summ += pow((data.F_xk_filtered[kol] - data.F_xk_noised[kol]) , 2);
											}
										temp_omega_delta.second = sqrt(summ/data.K_);

										temp_J = ( data.lambda[l] * temp_omega_delta.first) + (1 - data.lambda[l])*temp_omega_delta.second;

										if (temp_J < data.J)
											{
												data.J = temp_J;
												data.omega_delta =	temp_omega_delta;
												temp_save.alpha = temp_alpha;
												temp_save.metrics.delta = temp_omega_delta.second;
												temp_save.metrics.omega = temp_omega_delta.first;
												temp_save.metrics.J = temp_J;
												temp_save.graphic.Function_filtered = data.F_xk_filtered;
												temp_save.distance = sqrt(pow(temp_save.metrics.delta,2) + pow(temp_save.metrics.omega,2));
											}
									}
								data.J = 90000000.0;
								data.omega_delta = {90000000.0,90000000.0};
								save.push_back(temp_save);
							}
						
						if (jump == 0u)
							{
								v10.resut1.saves = save ;
								save.clear();
							}
						else
							{
								v10.resut2.saves = save ;
								save.clear();
							}
					}	
			
				find_best();
				std::cout<<"\t\t\t\t## For r = 3 ##\n";
				v10.resut1.print();
				std::cout<<"\t\t\t\t## For r = 5 ##\n";
				v10.resut2.print();
				/*std::cout<<"\n\n For r = 3\n";
				v10.resut1.best_save.graphic.print();
				std::cout<<"For r = 5 \n";
				v10.resut2.best_save.graphic.print();*/
			}
		void find_best()
			{
				Save temp_best;
				double best_distance = 999999900.0;
				for (Save save:v10.resut1.saves)
					{
						if (save.distance < best_distance)
							{
								best_distance = save.distance;
								temp_best = save;
							}
					}
				v10.resut1.best_save = temp_best;
				best_distance = 999999900.0;
				for (Save save:v10.resut2.saves)
					{
						if (save.distance < best_distance)
							{
								best_distance = save.distance;
								temp_best = save;
							}
					}
				v10.resut2.best_save = temp_best;
			}
	};
