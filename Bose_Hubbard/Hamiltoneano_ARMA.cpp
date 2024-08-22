#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <complex>
#include <armadillo>

using namespace std;

// Programamos la combinatoria C(n,k)=n!/k!(n-k)!
// recordemos que n>=k si k>n entonces C(n,k)=1
double combinatoria(int n, int k){
	int big; // Este va a ser el valor más grande entre k y n-k
	double n_des=n*1.0; // Aquí iremos desenciendo la n
	double result=1.0; // Aquí iremos calculando la combinatoria
	
	// Buscamos el valor más grande entre k y n-k
	if(n-k > k){
	big=(n-k);
	}else{
	big=k;
	}
	////////////
	
	for(int i = big; i >= 1; i--){ // En este for empezaremos a hacer esto (n/big)*(n-1/big-1)*...*(big)/(1)
		result=(n_des/i)*(result);
		n_des-=1;
	}
	return result;
}

// Nos busca un indice k tal que los numeros posteriores son 0
// k=3 si |2,1,0,0,0> bueno en c++ seria k=2
int searchk(int M,int n[]){
	int k=0; //iniciamos en ls posición k=0
	for(int i = 0; i<M; i++){ //Recorreremos todas las posiciones M para checar las j>i 
		k=i; //nuestro primer prospecto a ser k
		//cout<<"1 for "<<k<<endl;
		int flag=0; // definimos una bandera si es
			    // 1=> significa que encontro un n_{j} tal que los de adelante son 0
	                    // 0=> no lo ha encontrado
	                    
	        for(int j=k+1; j<M-1; j++){ //Revisaremos los n_{j} j>i si son cero
				//cout<<"2 for "<<j<<"="<<n[j]<<endl;
	        	if (n[j] != 0){ //si se encuentra una n_{j} != 0
	        		break; //ya no revises más
	        	}else if(j==M-2){ //si j llega al penultimo valor entonces 
	        		flag=1; // ya lo encontró
	        	}
	        	//cout<<"bandera "<<flag<<endl;
	        }	
	        if(flag==1){ //si la bandera es 1, es decir ya encontro la k
	        	k=i; //entonces esa i es la k
	        	break; //ya no hace falta revisar más
	        }
	        
	if (k==M-1){  //si no encontro ningún valor significa que tomaremos el penultimo valor
		k=M-2; // entonces k=M-1
	}
	}
	return k;
}

///Operadores del hamiltoneano///////////////////////////

//Operador ascenso, modifica vector
double ascenso_mu(int mu,int N,int vector_base[]){
	int part=vector_base[mu]+1;
	double norm=sqrt(vector_base[mu]+1);
	if(part>N){
		//part=0;
		norm=0.0;
	}
	vector_base[mu]=part;
	return norm;
}
//Operador descenso, modifica vector
double descenso_mu(int mu,int N,int vector_base[]){
	int part=vector_base[mu]-1;
	double norm=sqrt(vector_base[mu]);
	if(part<0){
		//part=0;
		norm=0.0;
	}
	vector_base[mu]=part;
	return norm;
}
//Operador numero, NO modifica vector
double Numero_mu(int mu,int N,int vector_base[]){
	return vector_base[mu];
}
////////////////////////////////////////////////////



//Funcion para localizar vectores base//////////
int Buscador(int vector_base[], int N, int M){


	int i0=N; // indice para ayudarnos a que este fijo el num de particulas disponibles
	int i1=0; //indice para ayudatnos a restar
	int Sitios=2;  //Numero de sitios que se quitan +1 (al principio se quita un sitio por eso es 2)
	int contador=0;
	
	
	for(int j=0; j<M; j++){
		//cout<<"INICIAMOS j="<<j<<endl;
		while(i0-i1 != vector_base[j]){	
		
			//cout<<"C("<<M-Sitios+i1<<","<<i1<<")="<<combinatoria(M-Sitios+i1,i1)<<endl;
			contador=contador+combinatoria(M-Sitios+i1,i1);
			i1++;
		}
		//cout<<"contador="<<contador<<endl;
		//cout<<"i1="<<i1<<endl;
		i0=i1;
		i1=0;
		Sitios++;
	}
	
	return contador;
	
	
}
//////////////////////////////////////////////////



//Funcion para localizar vectores base extendida//////////
int Buscador_base_ext(int vector_2_base_ext[], int N){	
	return vector_2_base_ext[0]*(N+1) + vector_2_base_ext[1];
}
//////////////////////////////////////////////////




int main(){
	int N=1; //Numero de particulas
	int M=3; //Numeros de sitios
	//N=2 y M=10
	//Parametros del Hamiltoniano
	double J=0.001;
	double K1=0;
	double K2=0;
	double U=1.0;
	
	//Lugares para trazar
	int n_nu=0;
	int n_mu=1;
	
	int D_exted=(N+1)*(N+1);
	//////////////////////////////////////////////////////////////
	//////CONTRUIMOS LA BASE DE FOCK//////////////////////////////
	/////////////////////////////////////////////////////////////
	int D=combinatoria(N+M-1,N); //declaramos la variable de la dimension del espacio de Hilbert
	//PARCHE BOBO
	double D_prueba=combinatoria(N+M-1,N);
	if (abs(D-D_prueba) < 0.0001){
		//cout<<"Holi_1"<<endl;
		D=combinatoria(N+M-1,N);
	}else{
		if (D > D_prueba){
			//cout<<"Holi_1.2"<<endl;
			D=D-1;
		}else if (D < D_prueba){
			//cout<<"Holi_1.2"<<endl;
			D=D+1;
		}
	}
	cout<<D<<" "<<D_exted<<endl;
	/////////////
	
	/*/Condicion inicial
	double Vector_0_F[D];
	//LLenamos la matriz de ceros
	for (int i=0; i<D; i++){
		Vector_0_F[i]=0;
	}
	/*/
	//int Vector_inicial[M]={0,0,1,0};
	//cout<<Buscador(Vector_inicial,N,M)<<endl;
	
	//Vector_0_F[0]=1.0; /// Vector de probas
	
	
	
	
	//Tiempo
	double dt=0.01;
	double Time_0=0;
	double Time_f=2500;
	double t;
	

	
	int DM[D][M]; //Aquí se colocaran las bases
	int n[M]; // Este vector va a ser el que modificaremos y usaremos para sobre escribir cada fila de DM
	// Llenamos la matriz de ceros
	for(int i = 0; i<D; i++){ 
		for(int j = 0; j<M; j++){
			DM[i][j]=0;
			n[j]=0;
		}
	}

	//cout<<"HOLA"<<endl;
	///////////////////
	///////////BASE EXTENDIDA/////////
	int base_2_extendida[D_exted][2];
	int l=0;
	for(int i = 0; i<N+1; i++){
		for(int j = 0; j<N+1; j++){
			base_2_extendida[l][0]=i;
			base_2_extendida[l][1]=j;
			l=l+1;
		}
	}
	//cout<<"HOLA"<<endl;
	//l=0;
	/*/ Imprimimos base extendida
	for(int i = 0; i<D_exted; i++){
		cout<<"["<<i<<"] ";
		for(int j = 0; j<2; j++){
			cout<< base_2_extendida[i][j]<<" \n"[j == 1];
		}
	}
	/*/
	//FIN DE BASE EXTENDIDA/////////
	/*/
	int zz[2]={3,3};
	cout<<Buscador_base_ext(zz,N)<<endl;
	/*/
	
	/////////////////////////////
	
	
	DM[0][0]=N; // El primer estado de la base
	n[0]=N; // Coincide con el estado base
	l=0; //NEcesitamos un contador para colocar los demas estados
	int k; // el indice donde despues de ese todos son ceros
	//cout<<serchk(M,n)<<endl;
	//cout<<"BASE DE FOCK"<<endl;
	while(n[M-1]!=N){
		l+=1; // incrementamos el contador
		k=searchk(M,n); //buscamos la k
		n[k]=n[k]-1; //Primer paso del cambio n_{k}=n_{k}-1
		//Compienza el segundo paso de la suma n_{k+1}=n_{k}-1
		n[k+1]=N;
		for(int i = 0; i<k+1; i++){ //se le coloca un +1 para que pueda acceder al menos a la posicion 0
			n[k+1]=n[k+1]-n[i];	
		}
		//Cambiamos por ceros los n_{i} tales que i >= k+2
		for(int i = k+2; i<M; i++){ 
			n[i]=0;	
		}
		//#Agregamos el estado n al la fila l
		for(int i=0; i<M; i++){
			DM[l][i]=n[i];
		}	
	}
	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////
	//////FIN DE LA CONSTRUCCION DE LA BASE DE FOCK//////////
	/*/ Imprimims la matriz
	for(int i = 0; i<D; i++){
		cout<<"["<<i<<"] ";
		for(int j = 0; j<M; j++){
			cout<< DM[i][j]<<" \n"[j == M-1];
		}
	}
	/*//////////////////
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	
	/////////////////////////////////////////////////////////////
	//////CONTRUIMOS LA MATRIZ HAMILTONEANA//////////////////////
	/////////////////////////////////////////////////////////////

	double a; // la constante que acompaña al operador ascenso 
	double at; // la constante que acompaña al operador descenso
	double num;// la constante que acompaña al operador numero
	int col;
	cout<<"HAMILTONIANO"<<endl;
	//cout<<D<<" "<<D_exted<<endl;
	//double H[D][D];
	cout<<D<<endl;
	arma::mat H(D,D);
	//cout<<"HAMILTONIANO"<<endl;
	/*/ Llenamos la matriz de ceros
	for(int i = 0; i<D; i++){ 
		for(int j = 0; j<D; j++){
			H[i][j]=0;
		}
	}
	/*/ 
	//AQui se empieza!!!!
	//cout<<"HAMILTONIANO"<<endl;
	
	if (J != 0){
		////// Esta es la parte     at_{m}*a_{m+1} + at_{m+1}*a_{m}
		for(int l1=0; l1<D; l1++){ //Usaremos cada elemento de la base
			for(int mu=0; mu<M; mu++){ // Recorreremos cada sitio del elemento de la base
				//construimos el vector con el que se trabajará
				for(int i = 0; i<M ; i++){
					n[i]=DM[l1][i];
				}
				//Aplicamos las condiciones de frontera periodicas (ANilli)
				if(mu == M-1){
					a=descenso_mu(0,N,n); // aplicamos el op de descenso en el comienzo
					at=ascenso_mu(mu,N,n); // aplicamos el op de ascenso
				}else{
					a=descenso_mu(mu+1,N,n); // aplicamos el op de descenso
					at=ascenso_mu(mu,N,n);  // aplicamos el op de ascenso
				}
				if (at*a != 0){
					col=Buscador(n, N, M);	
					//H[col][l1]+=(-J)*at*a;
					//H[l1][col]+=(-J)*at*a;
					H(col,l1)+=(-J)*at*a;
					H(l1,col)+=(-J)*at*a;
					
				}
			}
		}
	}if (K1 != 0){
		////// Esta es la parte     at_{mu}*n_{mu+nu}*a_{m+1} + at_{mu+1}*n_{mu+nu}*a_{m}
		for(int l1=0; l1<D; l1++){ //Usaremos cada elemento de la base
			for(int mu=0; mu<M; mu++){ // Recorreremos cada sitio del elemento de la base
				for(int nu=0; nu<2; nu++){ // esto es para el operador de numero n_{mu+nu}				
					//cout<<mu<<endl;
					//construimos el vector con el que se trabajará
					for(int i = 0; i<M ; i++){
						n[i]=DM[l1][i];
					}
					/// EStos IFS  SE PUEDE MEJORAR CON UNA FUNCION PERIODICA MIENTRAS SERÁ CON IFS
					if(mu == M-1 && nu == 1){
						a=descenso_mu(0,N,n); // aplicamos el op de descenso
						num=Numero_mu(0,N,n);// aplicamos el op de numero
						at=ascenso_mu(mu,N,n);  // aplicamos el op de ascenso	
					}else if(mu == M-1 && nu == 0){
						a=descenso_mu(0,N,n); // aplicamos el op de descenso
						num=Numero_mu(mu,N,n);// aplicamos el op de numero
						at=ascenso_mu(mu,N,n);  // aplicamos el op de ascenso
					}else{
						a=descenso_mu(mu+1,N,n); // aplicamos el op de descenso
						num=Numero_mu(mu+nu,N,n);// aplicamos el op de numero
						at=ascenso_mu(mu,N,n);  // aplicamos el op de ascenso
					}
					//cout<<at*num*a<<endl;
					if (at*num*a != 0){
						col=Buscador(n, N, M);
						//H[col][l1]+=(-K1)*at*num*a;
						//H[l1][col]+=(-K1)*at*num*a;
						H(col,l1)+=(-K1)*at*num*a;
						H(l1,col)+=(-K1)*at*num*a;
					}	
				}	
			}
		}
	}if (K2 != 0){
		////// Esta es la parte     at_{m}*at_{m}*a_{m+1}*a_{m+1} + at_{m+1}*at_{m+1}*a_{m}*a_{m}
		for(int l1=0; l1<D; l1++){ //Usaremos cada elemento de la base
			for(int mu=0; mu<M; mu++){ // Recorreremos cada sitio del elemento de la base
				//construimos el vector con el que se trabajará
				for(int i = 0; i<M ; i++){
					n[i]=DM[l1][i];
				}
				//Aplicamos las condiciones de frontera periodicas (ANilli)
				if(mu == M-1){
					a=descenso_mu(0,N,n); // aplicamos el op de descenso en el comienzo
					a=a*descenso_mu(0,N,n); // aplicamos el op de descenso en el comienzo
					
					at=ascenso_mu(mu,N,n); // aplicamos el op de ascenso
					at=at*ascenso_mu(mu,N,n); // aplicamos el op de ascenso
				}else{
					a=descenso_mu(mu+1,N,n); // aplicamos el op de descenso
					a=a*descenso_mu(mu+1,N,n); // aplicamos el op de descenso
					
					at=ascenso_mu(mu,N,n);  // aplicamos el op de ascenso
					at=at*ascenso_mu(mu,N,n);  // aplicamos el op de ascenso
				}
				if (at*a != 0){
					col=Buscador(n, N, M);	
					//H[col][l1]+=(-K2)*at*a;
					//H[l1][col]+=(-K2)*at*a;
					H(col,l1)+=(-K2)*at*a;
					H(l1,col)+=(-K2)*at*a;
					
				}
			}
		}
	}if (U != 0){
	////// Esta es la parte   n_{m}(n_{m}-1).. POTENCIAL
		for(int l1=0; l1<D; l1++){ //Usaremos cada elemento de la base
			//construimos el vector con el que se trabajará
			for(int i = 0; i<M ; i++){
				n[i]=DM[l1][i];
			}
			for(int i = 0; i<M ; i++){
				//H[l1][l1]+=U*0.5*(n[i]*(n[i]-1));
				H(l1,l1)+=U*0.5*(n[i]*(n[i]-1));
			}
			
		}	
	}	
	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////
	//////FIN DE LA CONSTRUCCION DE LA MATRIZ////////////////
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	/*/ Imprimims la matriz
	for(int i = 0; i<D; i++){
		cout<<"["<<i<<"] ";
		for(int j = 0; j<D; j++){
			cout<< H[i][j]<<"  "<<" \n"[j == D-1];
		}
	}
	/*/
	
	//cout<<"TERMINE EL HAMILTONIANO"<<endl;
	
	
	///////////////////////////////////////////////
	/////// APARTIR DE AQUI SE USA ARMADILLO///////
	///////////////////////////////////////////////
	//arma::mat H_A(&H[0][0], D, D); // PAsamos el H a armadillo H_A	
	arma::vec eigvals; // Aquí se colocan los eigenvalores
	arma::vec eigvals2;
	arma::mat Change_E_F; // Aqui se colocan los eigenvectores (columnas)
	eig_sym(eigvals,Change_E_F,H,"std"); // Se obtienen los eigen...
	arma::mat Change_F_E=Change_E_F.t(); // MAtriz cambio de base de Fock a Energias
	//arma::vec Vec0_F(&Vector_0_F[0],D); // Pasamos del vector inicial al vector armadillo
	arma::vec Vec0_F(D);
	Vec0_F(0)=1.0;
	arma::vec Vec0_E=Change_F_E*Vec0_F; //Cambiamos a la base de Energías
	arma::cx_vec Evol_temporal_E(D); // Declaramos el operador  evolución temporal
	arma::cx_vec Edo_F_t(D); //DEclaramos el estado al tiempo t en la base de energias
	arma::cx_mat Rho_C(D,D);
	arma::cx_mat Rho_C_reducida(D_exted,D_exted);
	arma::cx_mat Rho_C_trans_parcial(D_exted,D_exted);
	//cout<<sum(square(Vec0_E))<<endl;
	int bra;
	int ket;
	int NOSE[M];
	int base[M];
	int base_ex1[2];
	int base_ex2[2];	
	double Negativity;
	//H_A.print("Hamiltoniano");
	//cout.precision(17);
	//eigvals.raw_print(cout, "Eigenvalores");
	//eigvals.print("Eigenvalores");
	
	
	
	
	int Prueba=0;
	arma::vec Resultados(D);
	for(int i=0; i<D; i++){
		Resultados=(H*Change_E_F.col(i))-(eigvals(i)*Change_E_F.col(i));
		for(int j=0; j<D; j++){
			if (abs(Resultados(j))>0.000000000000001){
				Prueba+=1;
			}
		}
	}
	cout<<"Resultado de Prueba: "<<Prueba<<endl;
	
	
	
	
	
	
	///////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////AQUI EMPIEZA LA OTENCION DE LA NEGATIVIDAD///////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
	ofstream Matrix ("Negativity.dat"); //Al archivo que va a salir
	/*/
	for(int i = 0; i<D_exted; i++){
		cout<<"["<<i<<"] ";
		for(int j = 0; j<2; j++){
			cout<< base_2_extendida[i][j]<<" \n"[j == 1];
		}
	}
	/*/
	//cout<<"WHY|"<<base_2_extendida[0][0]<<","<<base_2_extendida[0][1]<<"><"<<base_2_extendida[1][0]<<","<<base_2_extendida[1][1]<<"|"<<endl;
	t=Time_0;
	//Vec0_F.print("Vector inicial en base de fock");
	while(t<Time_f){
		//cout<<"Entre al while al tiempo  "<<t<<endl;
		
		
		////////EVOLUCION TEMPORAL DEL ESTADO INICIAL////////////
		//t=0.5; // Tiempo
		Evol_temporal_E=exp((t)*(-1i)*eigvals);//Operador Evolución temporal 
		for(int i=0; i<D; i++){
			Edo_F_t(i)=Vec0_E(i)*Evol_temporal_E(i);
		}
		Edo_F_t=Change_E_F*Edo_F_t; //El cambio de base
		//cout<<"tiempo"<<t<<endl;
		//Edo_F_t.print("Estado evolucionado a tiempo");
		////////////////////////////////////////////////////////
		
		
		//////MATRIZ DE DENSIDAD////////////////////////////////
		Rho_C=Edo_F_t*Edo_F_t.t();
		//Rho_C.print("RHo normal");
		//cout<<trace(Rho_C)<<endl;
		///////////////////////////////////////////////////////
		
		

		
		//cout<<sum(square(Edo_F_t))<<endl;
		//H_A.print("Hamiltoniano");
		//eigvals.print("Eigenvalores");
		//Change_E_F.print("Eigenvectores");
		//Vec0_F.print("Vector inicial en base de fock");
		//Vec0_E.print("Vector inicial en base de Energias");
		//Edo_F_t.print("VEctoe evolucionado en base de fock");

		////MATRIZ DE DENSIDAD REDUCIDA EN LA BASE EXTENDIAD///////
		//Revisamos cada uno de los elemtnos de la base extendida |n_{\nu}^{i},n_{\mu}^{i}><n_{\nu}^{j},n_{\mu}^{j}|
		arma::cx_mat Rho_C_reducida(D_exted,D_exted);
		arma::cx_mat Rho_C_trans_parcial(D_exted,D_exted);
		for(int i=0; i<D_exted; i++){
			for(int j=0; j<D_exted; j++){
				//cout<<"Indices "<<i<<","<<j<<endl;
				//Esta es la primera delta de dirac (conservacion de particulas en la traza) --->
				//---> \delta_{ n_{\nu}^{i}+n_{\mu}^{i} , n_{\nu}^{j}+n_{\mu}^{j}
				bra=base_2_extendida[i][0]+base_2_extendida[i][1];
				ket=base_2_extendida[j][0]+base_2_extendida[j][1];
				//cout<<"1|"<<base_2_extendida[i][0]<<","<<base_2_extendida[i][1]<<"><"<<base_2_extendida[j][0]<<","<<base_2_extendida[j][1]<<"|"<<endl;
				if(bra==ket){
					//cout<<"HOLA__BRAKET "<<bra<<"="<<ket<<endl;
					//cout<<"|"<<base_2_extendida[i][0]<<","<<base_2_extendida[i][1]<<"><"<<base_2_extendida[j][0]<<","<<base_2_extendida[j][1]<<"|"<<endl;
					//Aqui es la segunda delta de dirac (Conservacion de particulas sin traza) --->
					//----->\delta_{n_{1}^{i}+n_{2}^{i}+...+n_{M}^{i}} , N}
					for(int k=0; k<D; k++){ //Tecnicamente ya sabemos quienes cumplen eso son los de la base de fock
					//cout<<"Estamos en "<<k<<endl;
					//reviamos que elementos de la base de fock tienen n_{\nu}^{i} y n_{\mu}^{i} simultaneamente de la base extendida para i
					//el k-esimo elemento es <k|\rho(t)| NOSE >
						if(DM[k][n_nu]==base_2_extendida[i][0] && DM[k][n_mu]==base_2_extendida[i][1]){
						//cout<<"Estamos en el if "<<k<<endl;
						//falta construir el NOSE que se compone de |n_{}^{i}...n_{\nu}^{j}...n_{\mu}^{j}>
		                		//Con el <k| coinciden los que no son n_{\nu}^{j} n_{\mu}^{j}, entonces
		                		for(int i1=0; i1<M; i1++){
		                			NOSE[i1]=DM[k][i1];
		                			base[i1]=DM[k][i1];
		                		}
		                		NOSE[n_nu]=base_2_extendida[j][0];
		                		NOSE[n_mu]=base_2_extendida[j][1];
		                		Rho_C_reducida(i,j)+=Rho_C(Buscador(base,N,M),Buscador(NOSE,N,M));
						}
					}
				}
			}
		}
		//cout<<"HOLAAAAAAAAAAAAAAA"<<endl;
		//Rho_C_reducida.print("Rho_Reducida");
		//Vamos a recorer cada elemento de la matriz de densidad reducida para transponerla
		for(int i=0; i<D_exted; i++){
			for(int j=0; j<D_exted; j++){
				//Reciclaremos variables bra será Transpo1
				//			 ket será Transpo2
				//base_ex[2]={base_2_extendida[i][0],base_2_extendida[j][1]};
				base_ex1[0]=base_2_extendida[i][0];
				base_ex1[1]=base_2_extendida[j][1];
				bra=Buscador_base_ext(base_ex1,N);
				//base_ex[2]={base_2_extendida[j][0],base_2_extendida[i][1]};
				base_ex2[0]=base_2_extendida[j][0];
				base_ex2[1]=base_2_extendida[i][1];
				ket=Buscador_base_ext(base_ex2,N);
				//cout<<i<<","<<j<<"--->"<<bra<<","<<ket<<endl;	
				Rho_C_trans_parcial(bra,ket)=Rho_C_reducida(i,j);
			}
		}
		//Rho_C_trans_parcial.print("Rho_reducida transpuesa parcial");
		//Rho_C_trans_parcial;
		eig_sym(eigvals2,Rho_C_trans_parcial);
		//eig_sym(eigvals2,Rho_C_reducida);
		//eigvals2.print("Eigenvalores_trans_partial");
		Negativity=0;
		for(int i=0; i<D_exted; i++){
			if(eigvals2(i)<0){
				Negativity+=abs(eigvals2(i));
			}
		}
		//cout<<"HOLA "<<Negativity<<endl;
		Matrix<<Negativity<<" "<<t<<endl;
		t=t+dt;
		
		//system("clear");
	}
	//////////////////
	Matrix.close();
	arma::arma_version ver;
	std::cout << "ARMA version: "<< ver.as_string() << std::endl;
	return 0;
	///////////////////
}
