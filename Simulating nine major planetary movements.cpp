#include<iostream>
#include<math.h>
#include<fstream>
#include<iomanip>
using namespace std;

double g = 6.67384e-11;     //���O�`�� 
double ms = 1.989e30;       //�Ӷ���q 
double dt;                  //�y��W�U�I���ɶ����� 

//�D�U�I��U�[�t�ת�x , y���q 
void fora(double& ax, double& ay)
{
	double x = ax;
	double y = ay;
	double r = pow(pow(x, 2)+pow(y, 2),0.5); //�U�I�P�Ӷ��Z�� 
	ax = -(g*ms/pow(r, 3))*x;   //	�U�I�V�ߥ[�t�ת�x , y���q  
	ay = -(g*ms/pow(r, 3))*y;
}

//�s�ծw��k 
//�N�t�׻P��m���O����x , y����q�å��s�ծw��k�p��U���I���� 
void RK4(double& x, double& y, double& vx, double& vy)
{
	double ax, ay;
	ax = x;
	ay = y;
	fora(ax, ay);           
	double vxk1 = dt*ax;   //�t�׵���ɶ����[�t�� 
	double vyk1 = dt*ay;
	double xk1 = dt*vx;    //�첾����ɶ����t�� 
	double yk1 = dt*vy;
	
	double xt = x+0.5*xk1;
	double yt = y+0.5*yk1;
	double vxt = vx+0.5*vxk1;
	double vyt = vy+0.5*vyk1;
	
	fora(xt, yt);		//���N�H�Wxt yt�ܴ���a
	double vxk2 = dt*xt;
	double vyk2 = dt*yt;
	double xk2 = dt*vxt;
	double yk2 = dt*vyt;
	
	xt = x+0.5*xk2;
	yt = y+0.5*yk2;
	vxt = vx+0.5*vxk2;
	vyt = vy+0.5*vyk2;
	
	fora(xt, yt);			//�Nxt yt�ഫa 
	double vxk3 = dt*xt;
	double vyk3 = dt*yt;
	double xk3 = dt*vxt;
	double yk3 = dt*vyt;
	
	xt = x+xk3;
	yt = y+yk3;
	vxt = vx+vxk3;
	vyt = vy+vyk3;

	fora(xt, yt);
	double vxk4 = dt*xt;
	double vyk4 = dt*yt;
	double xk4 = dt*vxt;
	double yk4 = dt*vyt;
	
	x = x + xk1/6.0 + xk2/3.0 + xk3/3.0 + xk4/6.0; 
	y = y + yk1/6.0 + yk2/3.0 + yk3/3.0 + yk4/6.0;
	vx = vx + vxk1/6.0 + vxk2/3.0 + vxk3/3.0 + vxk4/6.0;
	vy = vy + vyk1/6.0 + vyk2/3.0 + vyk3/3.0 + vyk4/6.0;
} 

//�ƾڿ�X�]�w 
void dataout(ofstream& fout, double t, double x, double y, double vx, double vy)
{
	fout << setw(15) << t 
		 << setw(15) << x 
		 << setw(15) << y 
		 << setw(15) << vx 
		 << setw(15) << vy << endl;
}

int main()
{
	cout << "�п�Jdt(s):" ;         //��J�ɶ����� , ��쬰�� 
	cin >> dt;
	double t, r, per, sum=0.0, ha, rd;
	double cross[6];
	// ���P
	ofstream a("Mercury.txt");
	a << setiosflags(ios::scientific);
	a << setw(15) << "#time"            // #sign is for gnuplot
	  << setw(15) << "x"
	  << setw(15) << "y"
	  << setw(15) << "vx" 
	  << setw(15) << "vy" << endl;
	double mervx=2.38e4, mervy=-3.29e4, merx=-4.79e10, mery=-4.57e10;//�]�w��l��m , �t��
	double min, max;//����Z�P���Z 
	double longr = 0.0 , dis = 0.0;
	t = 0.0;
	double y = mery;
	int i =0;
	max = pow(pow(merx, 2)+pow(mery, 2), 0.5);
	min = pow(pow(merx, 2)+pow(mery, 2), 0.5);
	//do�j��]�w��y��q�Ly�b , i++ , ��i�W�L6�h����{�� 
	do
	{
		t += dt;
		RK4(merx, mery, mervx, mervy);
		dataout(a, t, merx, mery, mervx, mervy);
		//cross[]�ΨӨD�g�� , ��y��q�Ly�b�h�O����U�ɶ� 
		if(y*mery<=0.0)
		{
			cross[i] = t;
			i++;
		}
		y = mery;
		r = pow(pow(merx, 2)+pow(mery, 2), 0.5);//�C���I�P�Ӷ��Z��
		if(r>=max)//�D����Z , �Q��if�j��D�C���I�P�Ӷ��Z�����̤j�� 
		{
			max = r; 
		}
		if(r<=min)//�D���Z , �Q��if�j��D�C���I�P�Ӷ��Z�����̤p�� 
		{
			min = r; 
		}
		
	
	} while(i<6);
	for(int n=0; n<5; n++)//
	{
		per = (cross[n+1]-cross[n])/86400/365*2;
		sum +=per;
	}
	cout << "���P����Z:" << max << endl;
	cout << "���P���Z:" << min << endl;
	cout << "���P¶��g��T(�~):" << sum/5 << endl; 
	ha = (max+min)/2/149597871000;
	cout << "���P¶��b���ba(AU):" << ha << endl;
	rd = pow(sum/5, 2)/pow(ha,3);
	cout << "�ĤT�w�� T^2/a^3 = " << rd << endl;  
	cout << "-----------------------------------------------------------------------------" << endl;
	a.close();
	
	//���P
	sum =0.0;
	ofstream b("Venus.txt");
	b << setiosflags(ios::scientific);
	b << setw(15) << "#time"            // #sign is for gnuplot
	  << setw(15) << "x"
	  << setw(15) << "y"
	  << setw(15) << "vx" 
	  << setw(15) << "vy" << endl;	
	double venvx=3.21e4, venvy=1.32e04, venx=4.18e10, veny=-1.00e11;
	t = 0.0;
	y = veny;
	i =0;
	min = pow(pow(venx, 2)+pow(veny, 2), 0.5);
	max = pow(pow(venx, 2)+pow(veny, 2), 0.5);
	do
	{
		t += dt;
		RK4(venx, veny, venvx, venvy);
		dataout(b, t, venx, veny, venvx, venvy);
		if(y*veny<=0.0)
		{
			cross[i] = t;
			i++;
		}
		y = veny;
		r = pow(pow(venx, 2)+pow(veny, 2), 0.5);
		if(r>=max)
		{
			max = r; 
		}
		if(r<=min)
		{
			min = r;
		}
	} while(i<6);
	for(int n=0; n<5; n++)
	{
		per = (cross[n+1]-cross[n])/86400/365*2;
		sum +=per;
	}
	cout << "���P����Z:" << max << endl;
	cout << "���P���Z:" << min << endl;
	cout << "���P¶��g��T(�~):" << sum/5 << endl; 
	ha = (max+min)/2/149597871000;
	cout << "���P¶��b���ba(AU):" << ha << endl;
	rd = pow(sum/5, 2)/pow(ha,3);
	cout << "�ĤT�w�� T^2/a^3 = " << rd << endl;  
	b.close();
	cout << "-----------------------------------------------------------------------------" << endl;
	
	//�a�y
	sum =0.0;
	ofstream c("Earth.txt");
	c << setiosflags(ios::scientific);
	c << setw(15) << "#time"            // #sign is for gnuplot
	  << setw(15) << "x"
	  << setw(15) << "y"
	  << setw(15) << "vx" 
	  << setw(15) << "vy" << endl;	
	double earvx=-3.01e04, earvy=3.13e03, earx=1.63e10, eary=1.47e11;
	t = 0.0;
	y = eary;
	i =0;
	min = pow(pow(earx, 2)+pow(eary, 2), 0.5);
	max = pow(pow(earx, 2)+pow(eary, 2), 0.5);
	do
	{
		t += dt;
		RK4(earx, eary, earvx, earvy);
		dataout(c, t, earx, eary, earvx, earvy);
		if(y*eary<=0.0)
		{
			cross[i] = t;
			i++;
		}
		y = eary;
		r = pow(pow(earx, 2)+pow(eary, 2), 0.5);
		if(r>=max)
		{
			max = r; 
		}
		if(r<=min)
		{
			min = r; 
		}
	} while(i<6);
	for(int n=0; n<5; n++)
	{
		per = (cross[n+1]-cross[n])/86400/365*2;
		sum +=per;
	}
	cout << "�a�y����Z:" << max << endl;
	cout << "�a�y���Z:" << min << endl;
	cout << "�a�y¶��g��T(�~):" << sum/5 << endl; 
	ha = (max+min)/2/149597871000;
	cout << "�a�y¶��b���ba(AU):" << ha << endl;
	rd = pow(sum/5, 2)/pow(ha,3);
	cout << "�ĤT�w�� T^2/a^3 = " << rd << endl;  
	c.close();
	cout << "-----------------------------------------------------------------------------" << endl;
	//���P
	sum =0.0;
	ofstream d("Mars.txt");
	d << setiosflags(ios::scientific);
	d << setw(15) << "#time"            // #sign is for gnuplot
	  << setw(15) << "x"
	  << setw(15) << "y"
	  << setw(15) << "vx" 
	  << setw(15) << "vy" << endl;	
	double marvx=2.10e4, marvy=-1.14e4, marx=-1.28e11, mary=-1.9e11;
	t = 0.0;
	y = mary;
	i =0;
	min = pow(pow(marx, 2)+pow(mary, 2), 0.5);
	max = pow(pow(marx, 2)+pow(mary, 2), 0.5);
	do
	{
		t += dt;
		RK4(marx, mary, marvx, marvy);
		dataout(d, t, marx, mary, marvx, marvy);
		if(y*mary<=0.0)
		{
			cross[i] = t;
			i++;
		}
		y = mary;
		r = pow(pow(marx, 2)+pow(mary, 2), 0.5);
		if(r>=max)
		{
			max = r; 
		}
		if(r<=min)
		{
			min = r; 
		}
	} while(i<6);
	for(int n=0; n<5; n++)
	{
		per = (cross[n+1]-cross[n])/86400/365*2;
		sum +=per;
	}
	cout << "���P����Z:" << max << endl;
	cout << "���P���Z:" << min << endl;
	cout << "���P¶��g��T(�~):" << sum/5 << endl; 
	ha = (max+min)/2/149597871000;
	cout << "���P¶��b���ba(AU):" << ha << endl;
	rd = pow(sum/5, 2)/pow(ha,3);
	cout << "�ĤT�w�� T^2/a^3 = " << rd << endl;  
	d.close();
	cout << "-----------------------------------------------------------------------------" << endl;
	//��P
	sum =0.0;
	ofstream e("Jupiter.txt");
	e << setiosflags(ios::scientific);
	e << setw(15) << "#time"            // #sign is for gnuplot
	  << setw(15) << "x"
	  << setw(15) << "y"
	  << setw(15) << "vx" 
	  << setw(15) << "vy" << endl;	
	double jupvx=1.13e4, jupvy=-5.66e3, jupx=-3.86e11, jupy=-7.04e11;
	t = 0.0;
	y = jupy;
	i =0;
	min = pow(pow(jupx, 2)+pow(jupy, 2), 0.5);
	max = pow(pow(jupx, 2)+pow(jupy, 2), 0.5);
	do
	{
		t += dt;
		RK4(jupx, jupy, jupvx, jupvy);
		dataout(e, t, jupx, jupy, jupvx, jupvy);
		if(y*jupy<=0.0)
		{
			cross[i] = t;
			i++;
		}
		y = jupy;
		r = pow(pow(jupx, 2)+pow(jupy, 2), 0.5);
		if(r>=max)
		{
			max = r; 
		}
		if(r<=min)
		{
			min = r; 
		}
	} while(i<6);
	for(int n=0; n<5; n++)
	{
		per = (cross[n+1]-cross[n])/86400/365*2;
		sum +=per;
	}
	cout << "��P����Z:" << max << endl;
	cout << "��P���Z:" << min << endl;
	cout << "��P¶��g��T(�~):" << sum/5 << endl; 
	ha = (max+min)/2/149597871000;
	cout << "��P¶��b���ba(AU):" << ha << endl;
	rd = pow(sum/5, 2)/pow(ha,3);
	cout << "�ĤT�w�� T^2/a^3 = " << rd << endl;  
	e.close();
	cout << "-----------------------------------------------------------------------------" << endl;
	//�g�P
	sum =0.0;
	ofstream f("Saturn.txt");
	f << setiosflags(ios::scientific);
	f << setw(15) << "#time"            // #sign is for gnuplot
	  << setw(15) << "x"
	  << setw(15) << "y"
	  << setw(15) << "vx" 
	  << setw(15) << "vy" << endl;	
	double satvx=-6.79e03, satvy=-7.39e03, satx=-1.04e12, saty=8.92e11;
	t = 0.0;
	y = saty;
	i =0;
	min = pow(pow(satx, 2)+pow(saty, 2), 0.5);
	max = pow(pow(satx, 2)+pow(saty, 2), 0.5);
	do
	{
		t += dt;
		RK4(satx, saty, satvx, satvy);
		dataout(f, t, satx, saty, satvx, satvy);
		if(y*saty<=0.0)
		{
			cross[i] = t;
			i++;
		}
		y = saty;
		r = pow(pow(satx, 2)+pow(saty, 2), 0.5);
		if(r>=max)
		{
			max = r; 
		}
		if(r<=min)
		{
			min = r; 
		}
	} while(i<6);
	for(int n=0; n<5; n++)
	{
		per = (cross[n+1]-cross[n])/86400/365*2;
		sum +=per;
	}
	cout << "�g�P����Z:" << max << endl;
	cout << "�g�P���Z:" << min << endl;
	cout << "�g�P¶��g��T(�~):" << sum/5 << endl; 
	ha = (max+min)/2/149597871000;
	cout << "�g�P¶��b���ba(AU):" << ha << endl;
	rd = pow(sum/5, 2)/pow(ha,3);
	cout << "�ĤT�w�� T^2/a^3 = " << rd << endl;  
	f.close();
	cout << "-----------------------------------------------------------------------------" << endl;
	//�Ѥ��P
	sum =0.0;
	ofstream g("Uranus.txt");
	g << setiosflags(ios::scientific);
	g << setw(15) << "#time"            // #sign is for gnuplot
	  << setw(15) << "x"
	  << setw(15) << "y"
	  << setw(15) << "vx" 
	  << setw(15) << "vy" << endl;	
	double uravx=1.86e3, uravy=6.22e3, urax=2.88e12, uray=-8.41e11;
	t = 0.0;
	y = uray;
	i =0;
	min = pow(pow(urax, 2)+pow(uray, 2), 0.5);
	max = pow(pow(urax, 2)+pow(uray, 2), 0.5);
	do
	{
		t += dt;
		RK4(urax, uray, uravx, uravy);
		dataout(g, t, urax, uray, uravx, uravy);
		if(y*uray<=0.0)
		{
			cross[i] = t;
			i++;
		}
		y = uray;
		r = pow(pow(urax, 2)+pow(uray, 2), 0.5);
		if(r>=max)
		{
			max = r; 
		}
		if(r<=min)
		{
			min = r; 
		}
	} while(i<6);
	for(int n=0; n<5; n++)
	{
		per = (cross[n+1]-cross[n])/86400/365*2;
		sum +=per;
	}
	cout << "�Ѥ��P����Z:" << max << endl;
	cout << "�Ѥ��P���Z:" << min << endl;
	cout << "�Ѥ��P¶��g��T(�~):" << sum/5 << endl; 
	ha = (max+min)/2/149597871000;
	cout << "�Ѥ��P¶��b���ba(AU):" << ha << endl;
	rd = pow(sum/5, 2)/pow(ha,3);
	cout << "�ĤT�w�� T^2/a^3 = " << rd << endl;  
	g.close();
	cout << "-----------------------------------------------------------------------------" << endl;
	//�����P 
	sum =0.0;
	ofstream h("Neptune.txt");
	h << setiosflags(ios::scientific);
	h << setw(15) << "#time"            // #sign is for gnuplot
	  << setw(15) << "x"
	  << setw(15) << "y"
	  << setw(15) << "vx" 
	  << setw(15) << "vy" << endl;	
	double nepvx=3.52e03, nepvy=4.14e3, nepx=3.4e12, nepy=-2.94e12;
	t = 0.0;
	y = nepy;
	i =0;
	min = pow(pow(nepx, 2)+pow(nepy, 2), 0.5);
	max = pow(pow(nepx, 2)+pow(nepy, 2), 0.5);
	do
	{
		t += dt;
		RK4(nepx, nepy, nepvx, nepvy);
		dataout(h, t, nepx, nepy, nepvx, nepvy);
		if(y*nepy<=0.0)
		{
			cross[i] = t;
			i++;
		}
		y = nepy;
		r = pow(pow(nepx, 2)+pow(nepy, 2), 0.5);
		if(r>=max)
		{
			max = r; 
		}
		if(r<=min)
		{
			min = r; 
		}
	} while(i<6);
	for(int n=0; n<5; n++)
	{
		per = (cross[n+1]-cross[n])/86400/365*2;
		sum +=per;
	}
	cout << "�����P����Z:" << max << endl;
	cout << "�����P���Z:" << min << endl;
	cout << "�����P¶��g��T(�~):" << sum/5 << endl; 
	ha = (max+min)/2/149597871000;
	cout << "�����P¶��b���ba(AU):" << ha << endl;
	rd = pow(sum/5, 2)/pow(ha,3);
	cout << "�ĤT�w�� T^2/a^3 = " << rd << endl;  
	h.close();
	
	system("pause");
}

