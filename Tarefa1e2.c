#include "mc920.h"

void getslice( char *filenameIN, char *c, int p, char *ocols, char *orows, char *q, char *filenameOUT ){
  	MedicalImage *I1 = NULL; /* Imagem de entrada */
  	GrayImage *outImage = NULL;
  	int  x, y;
  	int cols, rows;
  	rows = cols = 0;
	
	I1 = ReadMedicalImage( filenameIN ); //Le imagem medica -- funcao da biblioteca
		
	// condicoes para verificar qual a ordem de escrita	
	if( strcmp( c,"x" ) == 0 ){
		if( strcmp( ocols,"y" ) == 0 ){
			cols = I1->ny;
			rows = I1->nz;
		}
		else{
			cols = I1->nz;
			rows = I1->ny;		
		}		
	}
	else if( strcmp( c,"y" ) == 0 ){
		if( strcmp( ocols,"x" ) == 0 ){
			cols = I1->nx;
			rows = I1->nz;
		}
		else{
			cols = I1->nz;
			rows = I1->nx;		
		}	
	}
	else if( strcmp( c,"z" ) == 0 ){
		if( strcmp( ocols,"x" ) == 0 ){
			cols = I1->nx;
			rows = I1->ny;
		}
		else{
			cols = I1->ny;
			rows = I1->nx;		
		}		
	}
	else
		printf("\nErro..\n");
	
	
	outImage =  CreateGrayImage( cols, rows );	
	
	//direcao do corte for 00
	if( strcmp( q,"00" ) == 0 ){								
		for( y = 0; y < rows; y++){
			for( x = 0; x < cols; x++){
				if( strcmp( c,"x" ) == 0 ){ //verifica se eixo de corte é x
					if( strcmp( ocols,"y" ) == 0 ) //verifica a ordem y z ou z y
						outImage->val[y][x] = I1->val[y][x][p]; 
					else
						outImage->val[y][x] = I1->val[x][y][p]; 
				}					
				else if( strcmp( c,"y" ) == 0 ){ //verifica se eixo de corte é y
					if( strcmp( ocols,"x" ) == 0 ) //verifica a ordem x z ou z x
						outImage->val[y][x] = I1->val[y][p][x]; 
					else
						outImage->val[y][x] = I1->val[x][p][y];
				}
				else{ //se eixo de corte é z
					if( strcmp( ocols,"x" ) == 0 ) //verifica a ordem x y ou y x
						outImage->val[y][x] = I1->val[p][y][x]; //ra
					else
						outImage->val[y][x] = I1->val[p][x][y];
				}					
			}
		}
		WriteGrayImage( outImage, filenameOUT );			  			  				 
	}
	
	
	//direcao do corte for 01
	if( strcmp( q,"01" ) == 0 ){				
		for( y = 0; y < rows ; y++){
			for( x = 0; x < cols ; x++ ){
				if ( strcmp( c,"x" ) == 0 ){ //verifica se eixo de corte é x
					if( strcmp( ocols,"y" ) == 0 ) //verifica a ordem y z ou z y
						outImage->val[y][x] = I1->val[rows-1-y][x][p]; //ns
					else
						outImage->val[y][x] = I1->val[cols-1-x][y][p];
				}					
				else if( strcmp( c,"y" ) == 0 ){ //verifica se eixo de corte é y
					if( strcmp( ocols,"x" ) == 0 ) //verifica a ordem x z ou z x
						outImage->val[y][x] = I1->val[rows-1-y][p][x]; //nc-rc
					else
						outImage->val[y][x] = I1->val[cols-1-x][p][y];
				}
				else{ //se eixo de corte é z
					if( strcmp( ocols,"x" ) == 0 ) //verifica a ordem x y ou y x
						outImage->val[y][x] = I1->val[p][y][cols-1-x]; //na
					else
						outImage->val[y][x] = I1->val[p][x][rows-1-y];
				}
			}
		}
		WriteGrayImage( outImage, filenameOUT );			  				 
	}
	
	//direcao do corte for 10
	if( strcmp( q,"10" ) == 0 ){						
		for( y = 0; y < rows ; y++){
			for ( x = 0; x < cols; x++ ){      				
				if ( strcmp( c,"x" ) == 0 ){ //verifica se eixo de corte é x
					if ( strcmp( ocols,"y" ) == 0 ) //verifica a ordem y z ou z y
						outImage->val[y][x] = I1->val[y][cols-1-x][p];
					else
						outImage->val[y][x] = I1->val[x][rows-1-y][p];
				}						
				else if ( strcmp( c,"y" ) == 0 ){ //verifica se eixo de corte é y
					if ( strcmp( ocols,"x" ) == 0 ) //verifica a ordem x z ou z x
						outImage->val[y][x] = I1->val[y][p][cols-1-x];
					else
						outImage->val[y][x] = I1->val[x][p][rows-1-y];
				}
				else{ //se eixo de corte é z
					if ( strcmp( ocols,"x" ) == 0 ) //verifica a ordem x y ou y x
						outImage->val[y][x] = I1->val[p][rows-1-y][x];
					else
						outImage->val[y][x] = I1->val[p][cols-1-x][y];
				}
			}
		}
		WriteGrayImage( outImage, filenameOUT );			  			  				 
	}
	
	//direcao do corte for 11
	if( strcmp( q,"11" ) == 0 ){							    				
		for( y = 0; y < rows ; y++){
			for ( x = 0; x < cols; x++ ){
				if ( strcmp( c,"x" ) == 0 ){ //verifica se eixo de corte é x
					if ( strcmp( ocols,"y" ) == 0 ) //verifica a ordem y z ou z y
						outImage->val[y][x] = I1->val[rows-1-y][cols-1-x][p];
					else
						outImage->val[y][x] = I1->val[cols-1-x][rows-1-y][p];
				}					
				else if ( strcmp( c,"y" ) == 0 ){ //verifica se eixo de corte é y
					if ( strcmp( ocols,"x" ) == 0 ) //verifica a ordem x z ou z x
						outImage->val[y][x] = I1->val[rows-1-y][p][cols-1-x];
					else
						outImage->val[y][x] = I1->val[cols-1-x][p][rows-1-y];
				}
				else{ //se eixo de corte é z
					if ( strcmp( ocols,"x" ) == 0 ) //verifica a ordem x y ou y x
						outImage->val[y][x] = I1->val[p][rows-1-y][cols-1-x];
					else
						outImage->val[y][x] = I1->val[p][cols-1-x][rows-1-y];
				}
			}
		}
		WriteGrayImage( outImage, filenameOUT );		
  	}

}// fim funcao

void normaliza( char *filename, char *imageOUT ){
	int H;
	int  x, y, Imax, Imin;
	GrayImage *I;
	int aux;
	int nbits;
	

	I = ReadGrayImage( filename );
	
	Imax = INT_MIN; 
	Imin = INT_MAX;	
	for( y=0; y < I->ny; y++ ){
		for( x=0; x < I->nx; x++ ){
			if( I->val[y][x] > Imax )
				Imax = I->val[y][x];
		  	if( I->val[y][x] < Imin )
				Imin = I->val[y][x];
		}
	}
	
	nbits = 8;	
	H = pow( 2, nbits ) - 1;
	
	for( y=0; y < I->ny; y++ ){
		for( x=0; x < I->nx; x++ ){
			aux = I->val[y][x];
			if ( aux < Imin )
				aux = 0;
			else if( aux >= Imax )
				aux = H;
			else
				aux = (int) ( (H *(aux-Imin) ) / (Imax-Imin));
			I->val[y][x] = aux;	
		}
	}
	WriteGrayImage( I, imageOUT );
}


void negativo( char *filename, char *imageOUT ){
	int  x, y, Imax, Imin;
	GrayImage *I;
	int aux;
	
	
	I = ReadGrayImage( filename );
	
	Imax = INT_MIN; 
	Imin = INT_MAX;	
	for( y=0; y < I->ny; y++ ){
		for( x=0; x < I->nx; x++ ){
			if( I->val[y][x] > Imax )
				Imax = I->val[y][x];
		  	if( I->val[y][x] < Imin )
				Imin = I->val[y][x];
		}
	}
	
	for( y=0; y < I->ny; y++ ){
		for( x=0; x < I->nx; x++ ){
			aux = I->val[y][x];
			if ( aux < Imin )
				aux = Imax;
			else if( aux >= Imax )
				aux = Imin;
			else
				aux = (int) ( ((Imin-Imax) *(aux-Imin) / (Imax-Imin)) + Imax );	
			I->val[y][x] = aux;
		}
		
	}
	WriteGrayImage( I, imageOUT );
}


void bc( char *filename, float c, float b, char *imageOUT ){
	int H, nbits;
	int  x, y;
	int aux;
	float I1, I2;
	int Imax;
	float cont, bri;
	GrayImage *I;
	
	I = ReadGrayImage( filename );
	
	Imax = INT_MIN;	
	for( y=0; y < I->ny; y++ ){
		for( x=0; x < I->nx; x++ ){
			if( I->val[y][x] > Imax )
				Imax = I->val[y][x];
		}
	}
	
	if( Imax > 255 ){
		for( y = 0; y < I->ny; y++){
			for( x = 0; x < I->nx; x++){
				if( I->val[y][x] > Imax)
					I->val[y][x] = 255;                  
				else 
					I->val[y][x] = (int) ( ((float)I->val[y][x]/(float)Imax)*255 );                    
			}
		}
    }

	nbits = 8;		
	H = pow( 2, nbits ) - 1;
	
	//0-100
	cont = c/100;
	bri = b/100;
	
	//printf("B: %f\nC:%f", bri, cont);
	
	cont = (1-cont)*H;
	bri = (1-bri)*H; //inversão
	
	I1 = bri-(cont/2);
	I2 = bri+(cont/2);
	
	for( y=0; y < I->ny; y++ ){
		for( x=0; x < I->nx; x++ ){			
			aux = I->val[y][x];
			//printf("%d ", aux);
			if ( aux < I1 )
				aux = 0;
			else if( aux >= I2 )
				aux = 255;
			else
				aux = (int) ( H/(I2-I1) * (aux-I1) ); 
			I->val[y][x] = aux;				
			
			//printf("%d ", aux);
		}
	}
	WriteGrayImage( I, imageOUT );
}

void start();

void getsliceop(){
	int i;
	int op;
	int tam1, tam2;
	char filename[30], filenameOUT[30];
	char e[2], e1[2], e2[2];
	char q[3];
	
	char imageOUT[30];
	char imageOUTa[30];
	char imageOUTc[30];
	char imageOUTs[30];				
	
	do{
		printf("\nGetSlice.\n\n");
		printf("1. Visao do Radiologista.\n");
		printf("2. Visao do Neurocirurgiao.\n");
		printf("3. Outra visao.\n");
		printf("4. Voltar menu principal.\n");
		printf("0. Sair.\n");
		printf("\nDigite a opcao desejada: ");

		scanf("%d", &op);
		system("cls || clear");

		switch( op ){
			case 1:
				printf("Visao do Radiologista.\n\n");
		       	printf("Digite nome imagem de entrada(formato scn): ");
		          scanf("%s", filename );
		          tam1 = strlen( filename );
		          printf("Digite nome imagem de saida (sem formato): ");
		          scanf("%s", filenameOUT );
		          tam2 = strlen( filenameOUT );
		          printf("Digite posicao do corte: ");
		          scanf("%d", &i);
		          
		          filename[tam1] = '\0';
		          filenameOUT[tam2] = '\0';
		          
		          //Monta o nome do arquivo para salvar
				sprintf( imageOUTa, "%s-ra.pgm", filenameOUT );
				sprintf( imageOUTc, "%s-rc.pgm", filenameOUT );
				sprintf( imageOUTs, "%s-rs.pgm", filenameOUT );
		          
		          //printf("\n%s - %s\n", filename, imageOUTa );
				getslice( filename, "z", i, "x", "y", "00",  imageOUTa );
				getslice( filename, "y", i, "x", "z", "01",  imageOUTc );
				getslice( filename, "x", i, "y", "z", "11",  imageOUTs );
				
				system("cls || clear");
		          printf("Imagem Gerada com Sucesso!!!\n");
		          
				break;
            	case 2:
            		printf("Visao do Neurocirurgiao.\n\n");
		       	printf("Digite nome imagem de entrada(formato scn): ");
		          scanf("%s", filename );
		          tam1 = strlen( filename );
		          printf("Digite nome imagem de saida (sem formato): ");
		          scanf("%s", filenameOUT );
		          tam2 = strlen( filenameOUT );
		          printf("Digite posicao do corte: ");
		          scanf("%d", &i);
		          
		          filename[tam1] = '\0';
		          filenameOUT[tam2] = '\0';
		          
		          //Monta o nome do arquivo para salvar
				sprintf( imageOUTa, "%s-na.pgm", filenameOUT );
				sprintf( imageOUTc, "%s-nc.pgm", filenameOUT );
				sprintf( imageOUTs, "%s-ns.pgm", filenameOUT );
				
            		getslice( filename, "z", i, "x", "y", "01",  imageOUTa );
				getslice( filename, "y", i, "x", "z", "01",  imageOUTc );
				getslice( filename, "x", i, "y", "z", "01",  imageOUTs );
				
				system("cls || clear");
		          printf("Imagem Gerada com Sucesso!!!\n");
		          
		          break;
            	case 3:
            		printf("Outra visao.\n\n");
		       	printf("Digite nome imagem de entrada(formato scn): ");
		          scanf("%s", filename );
		          tam1 = strlen( filename );
		          printf("Digite nome imagem de saida (sem formato): ");
		          scanf("%s", filenameOUT );
		          tam2 = strlen( filenameOUT );
		          printf("Digite eixo do corte: ");
		          scanf("%s", e);
		          printf("Digite posicao do corte: ");
		          scanf("%d", &i);
		          printf("Digite proximo eixo: ");
		          scanf("%s", e1 );
		          printf("Digite proximo eixo: ");
		          scanf("%s", e2 );
		          printf("Digite quadrante de inicio:\nOpcoes: <00> <01> <10> <11>: ");
		          scanf("%s", q );
		          
		          filename[tam1] = '\0';
		          filenameOUT[tam2] = '\0';
		          
		          sprintf( imageOUT, "%s-%s.pgm", filenameOUT, q );
	
		          getslice( filename, e, i, e1, e2, q,  imageOUT );		          
		       
		          printf("Imagem Gerada com Sucesso!!!\n");
		          
                	break;  
                case 4:
                	start();
                	break;       
           	 case 0:
                	exit( -1 );
                	break;

            	default:
               	printf("Digite uma opcao valida.\n");
        	}
    	}while( op );
	
}

void start(){
	char filename[30], filenameOUT[30];;
	int op;
	int tam, tam2;
	float b, c;
	char imageOUT[30];
	
	do{		
		printf("\n1. GetSlice.\n");
		printf("2. Normalizar.\n");
		printf("3. Negativo.\n");
		printf("4. Brilho e Contraste.\n");
		printf("0. Sair\n");
		printf("\nDigite a opcao desejada: ");

		scanf("%d", &op);
		system("cls || clear");

		switch( op ){
			case 1:
				getsliceop();
				system("cls || clear");
		          printf("Imagem Gerada com Sucesso!!!\n");
		          
				break;
            	case 2:	
		       	printf("Normalizar.\n\n");
		       	printf("Digite nome imagem (formato pgm): ");
		          scanf("%s", filename );
		          tam = strlen( filename );
		          filename[tam] = '\0';
		          printf("Digite nome imagem de saida (sem formato): ");
		          scanf("%s", filenameOUT );
		          tam2 = strlen( filenameOUT );
		          filenameOUT[tam2] = '\0';          
		          
				sprintf( imageOUT, "%s-norm.pgm", filenameOUT );
		          normaliza( filename, imageOUT );
		          
		          system("cls || clear");
		          printf("Imagem Gerada com Sucesso!!!\n");
		          
		          filename[1]= '\0';
		          filenameOUT[1]= '\0';
		          break;
            	case 3:
               	printf("Negativo.\n\n");
		       	printf("Digite nome imagem (formato pgm): ");
		          scanf("%s", filename );
		          tam = strlen( filename );
		          filename[tam] = '\0';
		          printf("Digite nome imagem de saida (sem formato): ");
		          scanf("%s", filenameOUT );
		          tam2 = strlen( filenameOUT );
		          filenameOUT[tam2] = '\0';          
		          
				sprintf( imageOUT, "%s-neg.pgm", filenameOUT );						
		          negativo( filename, imageOUT );
		          
		          system("cls || clear");
		          printf("Imagem Gerada com Sucesso!!!\n");
		          
		          filename[1]= '\0';
		          filenameOUT[1]= '\0';
                	break;            
            	case 4:
            		printf("Brilho e Constraste.\n\n");
		       	printf("Digite nome imagem (formato pgm): ");
		          scanf("%s", filename );
		          tam = strlen( filename );
		          filename[tam] = '\0';
		          printf("Digite nome imagem de saida (sem formato): ");
		          scanf("%s", filenameOUT );
		          tam2 = strlen( filenameOUT );
		          filenameOUT[tam2] = '\0';
		          printf("Digite o brilho [0, 100]: ");
		          scanf("%f", &b);
		          printf("Digite constraste [0, 100]: ");
		          scanf("%f", &c);
		          
		          sprintf( imageOUT, "%s-b&c.pgm", filenameOUT );
		          
                	bc( filename, c, b, imageOUT );
                	
                	system("cls || clear");
		          printf("Imagem Gerada com Sucesso!!!\n");
		          
                	filename[1]= '\0';
		          filenameOUT[1]= '\0';
                	break;
           	 case 0:
                	exit( -1 );
                	break;

            	default:
               	printf("Digite uma opcao valida.\n");
        	}
    	}while( op );
}

int main(){
	
	start();
	
	return 0;
} //fim main

