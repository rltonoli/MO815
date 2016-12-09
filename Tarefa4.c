#include "mc920.h"

typedef struct voxel{
	float val[4]; // Real
}ValorVoxel;


// funcao sinal: se positivo retorna 1, se negativo retorna -1
int sign( int x ){
    if(x > 0)
    		return 1;
    if(x < 0) 
    		return -1;
    return 0;
}

//funcao desenhar reta entre dois pontos
void dda( float x1, float y1, float x2, float y2, char *filename, int H ){
	float p1[2];
	float p2[2];
	float p[2];
	int n, k;
	float du, dv;
	float Du, Dv;	
	GrayImage *img = NULL;
	
	FILE *fp;
	
	p1[0] = x1;
	p1[1] = y1;
	p2[0] = x2;
	p2[1] = y2;
	
	//printf("\nReta de p1 a p2:\n");
	//printf("x: %d y: %d - x: %d y: %d \n", p1[0], p1[1], p2[0], p2[1]);

	
	if( (p1[0] == p2[0]) && (p1[1] == p2[1]) ){
		n = 1;
	}
	else{
		Du = p2[0] - p1[0]; 
		Dv = p2[1] - p1[1];
		
		//printf("Du: %f Dv: %f\n", Du, Dv);
		
		if( abs(Du) >= abs(Dv) ){
			n = abs(Du)+1;
			du = sign(Du);
			dv = du * (Dv/Du);
			//printf("\ndv: %f\n", dv);
		}
		else{ 
			n = abs(Dv)+1;
			dv = sign(Dv);
			du = dv * (Du/Dv);
			//printf("\ndu: %f\n", du);;
		}
	}
	
	p[0] = p1[0];
	p[1] = p1[1];
	
	img = ReadGrayImage( filename );
	
	//printf("H: %d \n", H);

	fp = fopen(filename,"wb");
  
	for( k = 0; k < n; k++ ){		
		
		img->val[(int)(p[1])][(int)(p[0])] = H;
		
		//printf("p[0]: %d p[1]: %d \n", p[0], p[1]);
		p[0] = p[0] + du;
		p[1] = p[1] + dv;
				
	}
	
	WriteGrayImage( img, "wireframe.pgm");
	
	fclose(fp);
	
}


void rotacao( char *filename, int Rx, int Ry ){
     MedicalImage *I1 = NULL; /* Imagem de entrada */
     GrayImage *outImage = NULL; /* Imagem de saida*/
     int  x, y, z, i;
     int Imax;
     float tx, ty, tz;
     int D;
     float face;
     ValorVoxel faces[6];
     ValorVoxel result_fx[6];
     ValorVoxel result_fy[6];
	ValorVoxel v[8];
	ValorVoxel result_x[8];
	ValorVoxel result_y[8];
	
     I1 = ReadMedicalImage( filename ); //Le imagem medica -- funcao da biblioteca
     
     
     //Pegar valor maximo de intensidade
     Imax = 0;
     for (z=0; z < I1->nz; z++)
		for (y=0; y < I1->ny; y++)
		 	for (x=0; x < I1->nx; x++)
		 		if( I1->val[z][y][x] > Imax)
		 			Imax = I1->val[z][y][x];
     
     
     //printf("Imax: %d \n", I1->Imax);
     
     D = (int) sqrt( pow( I1->nx, 2 ) + pow( I1->ny, 2 ) + pow( I1->nz, 2 ) );
     
     
     
     tx = -(float)(I1->nx)/2;
     ty = -(float)(I1->ny)/2;
     tz = -(float)(I1->nz)/2;
     
     //printf("\nDados da imagem:\n");   
     //printf("Nx: %d Ny: %d Nz: %d\n", I1->nx, I1->ny, I1->nz );
     //printf("tx: %f ty: %f tz: %f\n\n", tx, ty, tz );
     //printf("D: %d\n", D);
     
     //translacao para o centro da coordenada dos 8 vertices
     v[0].val[0] = 0 + tx;
     v[0].val[1] = 0 + ty;
     v[0].val[2] = 0 + tz;
     v[0].val[3] = 0;
                           
     v[1].val[0] = I1->nx-1 + tx;
     v[1].val[1] = 0 + ty;
     v[1].val[2] = 0 + tz;
     v[1].val[3] = 0;
	              
     v[2].val[0] = 0 + tx;
     v[2].val[1] = I1->ny-1 + ty;
     v[2].val[2] = 0 + tz;
     v[2].val[3] = 0;

     v[3].val[0] = I1->nx-1 + tx;
     v[3].val[1] = I1->ny-1 + ty;
     v[3].val[2] = 0 + tz;
     v[3].val[3] = 0;
     
     v[4].val[0] = 0 + tx;
     v[4].val[1] = 0 + ty;
     v[4].val[2] = I1->nz-1 + tz;
     v[4].val[3] = 0;
                 
     v[5].val[0] = I1->nx-1+ tx;
     v[5].val[1] = 0 + ty;
     v[5].val[2] = I1->nz-1 + tz;
     v[5].val[3] = 0;
     
     v[6].val[0] = 0 + tx;
     v[6].val[1] = I1->ny-1 + ty;
     v[6].val[2] = I1->nz-1 + tz;
     v[6].val[3] = 0;	                    
                    
     v[7].val[0] = I1->nx-1 + tx;
     v[7].val[1] = I1->ny-1 + ty;
     v[7].val[2] = I1->nz-1 + tz;
     v[7].val[3] = 0;
	                
	
	/*
	printf("Translacao para centro da coordenada.\n");     
     for( i = 0; i < 8; i++ ){
          printf("i: %d >>> ", i);
          printf("x: %f ", v[i].val[0]);
          printf("y: %f ", v[i].val[1]);
          printf("z: %f ", v[i].val[2]);
          printf("C: %f\n", v[i].val[3]);
          
     }*/
     
     
     //printf( "\nCos: %f ", cos(Rx*PI/180) );
     //printf( "Sen: %f\n", sin(Rx*PI/180) );
     
        
     //rotacao em x 
     for( i = 0; i < 8; i++ ){      
          result_x[i].val[0] = v[i].val[0];
          result_x[i].val[1] = (v[i].val[1]*cos(Rx*PI/180)) - (v[i].val[2]*sin(Rx*PI/180));
          result_x[i].val[2] = (v[i].val[1]*sin(Rx*PI/180)) + (v[i].val[2]*cos(Rx*PI/180));         
     }
     
     /*
	printf("\nRotacao em Rx.\n");     
     for( i = 0; i < 8; i++ ){
          printf("i: %d >>> ", i);
          printf("x: %f ", result_x[i].val[0]);
          printf("y: %f ", result_x[i].val[1]);
          printf("z: %f ", result_x[i].val[2]);
          printf("C: %f\n", v[i].val[3]);         
     }*/
     
     //printf( "\nCos: %f ", cos(Ry*PI/180) );
     //printf( "Sen: %f\n", sin(Ry*PI/180) );
     
     //rotacao em y 
     for( i = 0; i < 8; i++ ){          
          result_y[i].val[0] = result_x[i].val[0]*cos(Ry*PI/180) + result_x[i].val[2]*sin(Ry*PI/180);
          result_y[i].val[1] = result_x[i].val[1];
          result_y[i].val[2] = -result_x[i].val[0]*sin(Ry*PI/180) + result_x[i].val[2]*cos(Ry*PI/180);
     } 
     
     /*
     printf("\nRotacao em Ry.\n");     
     for( i = 0; i < 8; i++ ){
          printf("i: %d >>> ", i);
          printf("x: %f ", result_y[i].val[0]);
          printf("y: %f ", result_y[i].val[1]);
          printf("z: %f ", result_y[i].val[2]);
          printf("C: %f\n", v[i].val[3]);         
     } */
     
     tx = D/2;
     ty = D/2;
     tz = D/2;

     
     //translacao para o centro da imagem
     for( i = 0; i < 8; i++ ){          
          v[i].val[0] = (int) (result_y[i].val[0] + tx);
          v[i].val[1] = (int) (result_y[i].val[1] + ty);
          v[i].val[2] = (int) (result_y[i].val[2] + tz);
     }
     
     /*
     printf("\nTranslacao para o centro da imagem:\n");      
     for( i = 0; i < 8; i++ ){
     	printf("i: %d >>> ", i);
          printf("x: %f ", v[i].val[0]);
          printf("y: %f ", v[i].val[1]);
          printf("z: %f ", v[i].val[2]);
          printf("C: %f\n", v[i].val[3]);
     }*/
     
     //atribuindo valores as faces
     i = 0;       
     if( i == 0 ){
	     faces[i].val[0] = 1;
	     faces[i].val[1] = 0;
	     faces[i].val[2] = 0;
	     faces[i].val[3] = 0;
	     i++;
	}
	if( i == 1 ){
		faces[i].val[0] = -1;
	     faces[i].val[1] = 0;
	     faces[i].val[2] = 0;
	     faces[i].val[3] = 0;
	     i++;		
	}
	if( i == 2 ){
	     faces[i].val[0] = 0;
	     faces[i].val[1] = 1;
	     faces[i].val[2] = 0;
	     faces[i].val[3] = 0;
	     i++;
	}
	if( i == 3 ){
		faces[i].val[0] = 0;
	     faces[i].val[1] = -1;
	     faces[i].val[2] = 0;
	     faces[i].val[3] = 0;
	     i++;		
	}
	if( i == 4 ){
	     faces[i].val[0] = 0;
	     faces[i].val[1] = 0;
	     faces[i].val[2] = 1;
	     faces[i].val[3] = 0;
	     i++;
	}
	if( i == 5 ){
		faces[i].val[0] = 0;
	     faces[i].val[1] = 0;
	     faces[i].val[2] = -1;
	     faces[i].val[3] = 0;
	     i++;		
	}
      
          
     //rotacao em x das faces
     for( i = 0; i < 6; i++ ){          
          result_fx[i].val[0] = faces[i].val[0];
          result_fx[i].val[1] = faces[i].val[1]*cos(Rx*PI/180) - faces[i].val[2]*sin(Rx*PI/180);
          result_fx[i].val[2] = faces[i].val[1]*sin(Rx*PI/180) + faces[i].val[2]*cos(Rx*PI/180);
     } 
     
     /*
     printf("\nRotacao em X face:\n");      
     for( i = 0; i < 6; i++ ){
     	printf("i: %d >>> ", i);
          printf("x: %f ", result_fx[i].val[0]);
          printf("y: %f ", result_fx[i].val[1]);
          printf("z: %f \n", result_fx[i].val[2]);
     }*/
     
     
     //rotacao em y das faces
     for( i = 0; i < 6; i++ ){          
          result_fy[i].val[0] = result_fx[i].val[0]*cos(Ry*PI/180) + result_fx[i].val[2]*sin(Ry*PI/180);
          result_fy[i].val[1] = result_fx[i].val[1];
          result_fy[i].val[2] = -result_fx[i].val[0]*sin(Ry*PI/180) + result_fx[i].val[2]*cos(Ry*PI/180);
     }
     
     /*
     printf("\nRotacao em Y face:\n");      
     for( i = 0; i < 6; i++ ){
     	printf("i: %d >>> ", i);
          printf("x: %f ", result_fy[i].val[0]);
          printf("y: %f ", result_fy[i].val[1]);
          printf("z: %f \n", result_fy[i].val[2]);
     }*/
     
     outImage =  CreateGrayImage( D, D );
          
     for( y = 0; y < D; y++ ){
     	for( x = 0; x < D; x++ ) 
     		outImage->val[y][x] = 0; // pintando imagem de preta	
	}
	
	WriteGrayImage( outImage, "wireframe.pgm");
     
     //produto interno 
     for( i = 0; i < 6; i++ ){          
          face =  0 * result_fy[i].val[0];
          face = face + ( 0 * result_fy[i].val[1] );
          face = face + ( -1 * result_fy[i].val[2] );
          
          //printf("\nFace: %f", face);
         
          if( face > 0){		     
		     
		     printf("\nFace visivel: ");
		     /*
		     printf("%f ",result_fy[i].val[0]);
		     printf("%f ",result_fy[i].val[1]);
		     printf("%f ",result_fy[i].val[2]);
		     printf("%f \n",v[i].val[3]);	*/	   
		     
		     // (-1,0,0,0) ou (1,0,0,0) 
		     if( i == 0 || i == 1 ){
		     	if( i == 1 ){ // os 4 pontos de cada face
		     		
		     		printf("-1 ");
					printf("0 ");
					printf("0 ");
					printf("0\n");
					
		     		dda( v[0].val[0], v[0].val[1], v[4].val[0], v[4].val[1], "wireframe.pgm", Imax );
		     		
		     		dda( v[4].val[0], v[4].val[1], v[6].val[0], v[6].val[1], "wireframe.pgm", Imax );	
		     		
		     		dda( v[6].val[0], v[6].val[1], v[2].val[0], v[2].val[1], "wireframe.pgm", Imax );
				
   					dda( v[2].val[0], v[2].val[1], v[0].val[0], v[0].val[1], "wireframe.pgm", Imax );
		     		
				}
				else{ // os 4 pontos de cada face
					
					printf("1 ");
					printf("0 ");
					printf("0 ");
					printf("0\n");
					
					dda( v[1].val[0], v[1].val[1], v[5].val[0], v[5].val[1], "wireframe.pgm", Imax );
		     		
		     		dda( v[5].val[0], v[5].val[1], v[7].val[0], v[7].val[1], "wireframe.pgm", Imax );	
		     		
		     		dda( v[7].val[0], v[7].val[1], v[3].val[0], v[3].val[1], "wireframe.pgm", Imax );
				
   					dda( v[3].val[0], v[3].val[1], v[1].val[0], v[1].val[1], "wireframe.pgm", Imax );
				}
			}			
			
			// (0,-1,0,0) ou (0,1,0,0) 
		     if( i == 2 || i == 3 ){
		     	if( i == 3 ){ // os 4 pontos de cada face
		     		
		     		printf("0 ");
					printf("-1 ");
					printf("0 ");
					printf("0\n");
		     		
		     		dda( v[0].val[0], v[0].val[1], v[1].val[0], v[1].val[1], "wireframe.pgm", Imax );
		     		
		     		dda( v[1].val[0], v[1].val[1], v[5].val[0], v[5].val[1], "wireframe.pgm", Imax );	
		     		
		     		dda( v[5].val[0], v[5].val[1], v[4].val[0], v[4].val[1], "wireframe.pgm", Imax );
				
   					dda( v[4].val[0], v[4].val[1], v[0].val[0], v[0].val[1], "wireframe.pgm", Imax );
						
				}
				else{ // os 4 pontos de cada face
					
					printf("0 ");
					printf("1 ");
					printf("0 ");
					printf("0\n");
					
					dda( v[2].val[0], v[2].val[1], v[3].val[0], v[3].val[1], "wireframe.pgm", Imax );
		     		
		     		dda( v[3].val[0], v[3].val[1], v[7].val[0], v[7].val[1], "wireframe.pgm", Imax );	
		     		
		     		dda( v[7].val[0], v[7].val[1], v[6].val[0], v[6].val[1], "wireframe.pgm", Imax );
				
   					dda( v[6].val[0], v[6].val[1], v[2].val[0], v[2].val[1], "wireframe.pgm", Imax );
				}
			}
						
		     // (0,0,-1,0) ou (0,0,1,0) 
		     if( i == 4 || i == 5 ){
		     	if( i == 5 ){
		     		printf("0 ");
					printf("0 ");
					printf("-1 ");
					printf("0\n");		     		
   		     		
		     		dda( v[0].val[0], v[0].val[1], v[1].val[0], v[1].val[1], "wireframe.pgm", Imax );
		     		
		     		dda( v[1].val[0], v[1].val[1], v[3].val[0], v[3].val[1], "wireframe.pgm", Imax );		     		
		     		
		     		dda( v[3].val[0], v[3].val[1], v[2].val[0], v[2].val[1], "wireframe.pgm", Imax );
		     		
		     		dda( v[2].val[0], v[2].val[1], v[0].val[0], v[0].val[1], "wireframe.pgm", Imax );
					
				}
				else{
					printf("0 ");
					printf("0 ");
					printf("1 ");
					printf("0\n");
					
					
					dda( v[4].val[0], v[4].val[1], v[5].val[0], v[5].val[1], "wireframe.pgm", Imax );					
					
   					dda( v[5].val[0], v[5].val[1], v[7].val[0], v[7].val[1], "wireframe.pgm", Imax ); 		     		
		     		
		     		dda( v[7].val[0], v[7].val[1], v[6].val[0], v[6].val[1], "wireframe.pgm", Imax );
		     		
		     		dda( v[6].val[0], v[6].val[1], v[4].val[0], v[4].val[1], "wireframe.pgm", Imax ); 
				}
			}
					     		     
		}
     }              	          
	          
}


int main(int argc, char **argv){	 
	
	int rx, ry;
	if( argc != 3 ){
		fprintf( stderr, "use: rotacao\n");
		exit( -1 );
	}
	
	rx = atoi( argv[1] );
	ry = atoi( argv[2] );
	
	
     rotacao( "skull.scn", rx, ry ); // dando erro com 45, 45...
	return 0;
} //fim main


