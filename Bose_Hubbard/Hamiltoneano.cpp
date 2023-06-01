#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
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




int main(){
	int N=2; //Numero de particulas
	int M=50; //Numeros de sitios
	
	double J=0;
	double K1=0;
	double K2=0;
	double U=1;
	//////////////////////////////////////////////////////////////
	//////CONTRUIMOS LA BASE DE FOCK//////////////////////////////
	/////////////////////////////////////////////////////////////
	int D=combinatoria(N+M-1,N); //declaramos la variable de la dimension del espacio de Hilbert
	//PARCHE BOBO
	double D_prueba=combinatoria(N+M-1,N);
	if (abs(D-D_prueba) < 0.00001){
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
	/////////////
	
	//cout<<D<<endl;
	//cout<<"::::::::::::::::::"<<endl;
	int DM[D][M]; //Aquí se colocaran las bases
	int n[M]; // Este vector va a ser el que modificaremos y usaremos para sobre escribir cada fila de DM
	
	
	// Llenamos la matriz de ceros
	for(int i = 0; i<D; i++){ 
		for(int j = 0; j<M; j++){
			DM[i][j]=0;
			n[j]=0;
		}
	}
	/////////////////////////////
	
	
	DM[0][0]=N; // El primer estado de la base
	n[0]=N; // Coincide con el estado base
	
	int l=0; //NEcesitamos un contador para colocar los demas estados
	int k; // el indice donde despues de ese todos son ceros
	//cout<<serchk(M,n)<<endl;
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

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	
	/////////////////////////////////////////////////////////////
	//////CONTRUIMOS LA MATRIZ HAMILTONEANA//////////////////////
	/////////////////////////////////////////////////////////////
	
	
	
	
	
	//l=37; // Reesiclamos y esta l lo usaremos para que num de vector de la base
	//k=2; // El sitio 
	
	double a; // la constante que acompaña al operador ascenso 
	double at; // la constante que acompaña al operador descenso
	double num;// la constante que acompaña al operador numero
	
	int col;
	cout<<"HAMILTONIANO_antes"<<endl;
	double H[D][D];
	cout<<"HAMILTONIANO_despues"<<endl;	
	// Llenamos la matriz de ceros
	for(int i = 0; i<D; i++){ 
		for(int j = 0; j<D; j++){
			H[i][j]=0;
		}
	}
	
	// AQui se empieza!!!!

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
				H[col][l1]+=(-J)*at*a;
				H[l1][col]+=(-J)*at*a;
				
			}
		}
	}

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
					H[col][l1]+=(-K1)*at*num*a;
					H[l1][col]+=(-K1)*at*num*a;
				}	
			}	
		}
	}
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
				H[col][l1]+=(-K2)*at*a;
				H[l1][col]+=(-K2)*at*a;
				
			}
		}
	}
	
	////// Esta es la parte   n_{m}(n_{m}-1)
	for(int l1=0; l1<D; l1++){ //Usaremos cada elemento de la base
		//construimos el vector con el que se trabajará
		for(int i = 0; i<M ; i++){
			n[i]=DM[l1][i];
		}
		for(int i = 0; i<M ; i++){
			H[l1][l1]+=U*0.5*(n[i]*(n[i]-1));
		}
		
	}
	
	
	// Imprimims la matriz
	for(int i = 0; i<D; i++){
		cout<<"["<<i<<"] ";
		for(int j = 0; j<D; j++){
			cout<< H[i][j]<<"  "<<" \n"[j == D-1];
		}
	}
	//////////////////
	// Salvar la matriz
	/*/
	ofstream Matrix ("Hamilt.dat"); //Al archivo que va a salir
	for(int i = 0; i<D; i++){
		for(int j = 0; j<D; j++){
			Matrix<<H[i][j]<<" \n"[j == D-1];
		}
	}
	Matrix.close();
	/*/
	return 0;
	
	
	
	///////////////////
}
