#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
using namespace std;

// Programamos la combinatoria C(n,k)=n!/k!(n-k)!
// recordemos que n>=k si k>n entonces C(n,k)=1
float combinatoria(int n, int k){
	int big; // Este va a ser el valor más grande entre k y n-k
	float n_des=n*1.0; // Aquí iremos desenciendo la n
	float result=1.0; // Aquí iremos calculando la combinatoria
	
	// Buscamos el valor más grande entre k y n-k
	if(n-k > k){
	big=(n-k);
	}else{
	big=k;
	}
	////////////
	
	for(int i = big; i >= 1; i--){ // En este for empezaremos a hacer esto (n/big)*(n-1/big-1)*...*(big)/(1)
		
		result=(n_des/i)*(result);
		
		n_des=n_des-1;
		//cout<<result<<endl;
	}
	//return ceil(result); //funcion techo para arreglar lo de float to int
	return result;
}

// Nos busca un indice k tal que los numeros posteriores son 0
// k=3 si |2,1,,0,0,0> bueno en c++ seria k=2
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


// Contruimos la base de fock
//N--> numero de particulas y M--> numero de sitios
int base_fock(int N, int M){ 
	int D=combinatoria(N+M-1,N); //declaramos la variable de la dimension del espacio de Hilbert
	float D_prueba=combinatoria(N+M-1,N);
	
	//cout<<D<<endl;
	//cout<<D_prueba<<endl;
	//cout<<D-D_prueba<<endl;
	//cout<<abs(D-D_prueba)<<endl;
	//cout<<0.00001<<endl;
	//PARCHE BOBO
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
	
	
	if (D < combinatoria(N+M-1,N)){
		cout<<"Holi"<<endl;
	}
	//cout<<D<<endl;
	
	int DM[D][M]; //Aquí se colocaran las bases
	int n[M]; // Este vector va a ser el que modificaremos y usaremos para sobre escribir cada fila de DM
	
	// Llenamos la matriz de ceros
	for(int i = 0; i<D; i++){ 
		for(int j = 0; j<M; j++){
			DM[i][j]=0;
			n[j]=0;
		}
	}
	///////////
	
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
	
	
	
	
	
	// Imprimims la matriz
	for(int i = 0; i<D; i++){
		cout<<"["<<i<<"] ";
		for(int j = 0; j<M; j++){
			cout<< DM[i][j]<<" \n"[j == M-1];
		}
	}
	///////////////////
	
	
	return 0;
}
/*
int index_periodica(int index, int M){
	if (index > M-1){
		index=0+abs(index-M-1);
	}
}
*/

int main(){
	base_fock(2,10);
	//base_fock(5,3); //Problemas con uno de menos celi
	//base_fock(3,7); //Problemas con tener uno de más
	//base_fock(3,7); //Problemas // N,M
	//cout<<combinatoria(5+3-1,5)<<endl;
	//cout<<combinatoria(5+3-1,5)<<endl;
	//cout<<D<<endl;
	
	
	//int D=combinatoria(5+3-1,5); //declaramos la variable de la dimension del espacio de Hilbert
	//float D_prueba=combinatoria(5+3-1,5);
	//cout<<D<<" "<<D_prueba<<" Diferencia "<<D_prueba-D<<endl;
	
	return 0;
}
