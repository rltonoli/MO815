#include "mainwindow.h"
#include <QApplication>
#include "mc920.h"

void getslice( char *filenameIN, char *c, int p, char *ocols, char *orows, char *q, char *filenameOUT ){
    MedicalImage *I1 = NULL; /* Imagem de entrada */
    GrayImage *outImage = NULL;
    int  x, y, z, i, j;
    int cols, rows;
    int Imax=0;
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
    else{ printf("\nErro..\n");
    }

    outImage =  CreateGrayImage( cols, rows );

    //quando o quadrante do corte for 00
    if( strcmp( q,"00" ) == 0 ){
        for ( y = 0; y < rows; y++){
            for ( x = 0; x < cols; x++){
                if ( strcmp( c,"x" ) == 0 ){ //verifica se eixo de corte é x
                    if ( strcmp( ocols,"y" ) == 0 ) //verifica a ordem y z ou z y
                        outImage->val[y][x] = I1->val[y][x][p];
                    else
                        outImage->val[y][x] = I1->val[x][y][p];
                }

                else if ( strcmp( c,"y" ) == 0 ){ //verifica se eixo de corte é y
                    if ( strcmp( ocols,"x" ) == 0 ) //verifica a ordem x z ou z x
                        outImage->val[y][x] = I1->val[y][p][x];
                    else
                        outImage->val[y][x] = I1->val[x][p][y];
                }
                else{ //se eixo de corte é z
                    if ( strcmp( ocols,"x" ) == 0 ) //verifica a ordem x y ou y x
                        outImage->val[y][x] = I1->val[p][y][x];
                    else
                        outImage->val[y][x] = I1->val[p][x][y];
                }

            }
        }
        WriteGrayImage( outImage, filenameOUT );
    }
    //quando o quadrante do corte for 01
    if( strcmp( q,"01" ) == 0 ){

        for ( y = 0; y < rows ; y++)
            for ( x = cols-1, i = 0; x > 0; x--, i++){
                if ( strcmp( c,"x" ) == 0 ){ //verifica se eixo de corte é x
                    if ( strcmp( ocols,"y" ) == 0 ) //verifica a ordem y z ou z y
                        outImage->val[y][x] = I1->val[y][x][p];
                    else
                        outImage->val[y][x] = I1->val[x][y][p];
                }

                else if ( strcmp( c,"y" ) == 0 ){ //verifica se eixo de corte é y
                    if ( strcmp( ocols,"x" ) == 0 ) //verifica a ordem x z ou z x
                        outImage->val[y][x] = I1->val[y][p][x];
                    else
                        outImage->val[y][x] = I1->val[x][p][y];
                }
                else{ //se eixo de corte é z
                    if ( strcmp( ocols,"x" ) == 0 ) //verifica a ordem x y ou y x
                        outImage->val[y][x] = I1->val[p][y][x];
                    else
                        outImage->val[y][x] = I1->val[p][x][y];
                }
            }
        WriteGrayImage( outImage, filenameOUT );
    }
    //quando o quadrante do corte for 10
    if( strcmp( q,"10" ) == 0 ){

        for ( y = rows-1, i = 0; y >= 0; y--, i++ )
            for ( x = 0; x < cols; x++ ){
                if ( strcmp( c,"x" ) == 0 ){ //verifica se eixo de corte é x
                    if ( strcmp( ocols,"y" ) == 0 ) //verifica a ordem y z ou z y
                        outImage->val[y][x] = I1->val[y][x][p];
                    else
                        outImage->val[y][x] = I1->val[x][y][p];
                }

                else if ( strcmp( c,"y" ) == 0 ){ //verifica se eixo de corte é y
                    if ( strcmp( ocols,"x" ) == 0 ) //verifica a ordem x z ou z x
                        outImage->val[y][x] = I1->val[y][p][x];
                    else
                        outImage->val[y][x] = I1->val[x][p][y];
                }
                else{ //se eixo de corte é z
                    if ( strcmp( ocols,"x" ) == 0 ) //verifica a ordem x y ou y x
                        outImage->val[y][x] = I1->val[p][y][x];
                    else
                        outImage->val[y][x] = I1->val[p][x][y];
                }
            }
        WriteGrayImage( outImage, filenameOUT );
    }
    //quando o quadrante do corte for 11
    if( strcmp( q,"11" ) == 0 ){
        for ( y = rows-1, i = 0; y > 0 ; y--, i++)
            for ( x = cols-1, j = 0; x > 0; x--, j++){
                if ( strcmp( c,"x" ) == 0 ){ //verifica se eixo de corte é x
                    if ( strcmp( ocols,"y" ) == 0 ) //verifica a ordem y z ou z y
                        outImage->val[y][x] = I1->val[y][x][p];
                    else
                        outImage->val[y][x] = I1->val[x][y][p];
                }

                else if ( strcmp( c,"y" ) == 0 ){ //verifica se eixo de corte é y
                    if ( strcmp( ocols,"x" ) == 0 ) //verifica a ordem x z ou z x
                        outImage->val[y][x] = I1->val[y][p][x];
                    else
                        outImage->val[y][x] = I1->val[x][p][y];
                }
                else{ //se eixo de corte é z
                    if ( strcmp( ocols,"x" ) == 0 ) //verifica a ordem x y ou y x
                        outImage->val[y][x] = I1->val[p][y][x];
                    else
                        outImage->val[y][x] = I1->val[p][x][y];
                }
            }
        WriteGrayImage( outImage, filenameOUT );
    }
}


float Maximo( float num1, float num2 ) {
    if( num1 > num2 )
        return num1;
    return num2;
}

void composicaoImagem( GrayImage *iGray, ColorImage *cI ){
    float halfH;
    int channel, x, y, r, g, b, Yi, Cg, Co, Imax;
    ColorImage *iAux = NULL, *iFinal = NULL;

    iAux = CreateColorImage( cI->nx, cI->ny );
    iFinal = CreateColorImage( cI->nx, cI->ny );

    Imax = INT_MIN;
    for( y=0; y < cI->ny; y++ ){
        for( x=0; x < cI->nx; x++ ){
            for( channel = 0; channel < 3; channel++ )
                if( cI->cor[y][x].val[channel] > Imax )
                    Imax = cI->cor[y][x].val[channel];
        }
    }


    // Transforma para CgCo recebendo o Y do slice da imagem médica
    halfH = Imax / 2;
    for (y = 0; y < cI->ny; y++) {
        for (x = 0; x < cI->nx; x++) {
            r = cI->cor[y][x].val[2];
            g = cI->cor[y][x].val[1];
            b = cI->cor[y][x].val[0];


            if( iGray->val[y][x] != 0 )
                iAux->cor[y][x].val[2] = iGray->val[y][x]; // Y

            else
                iAux->cor[y][x].val[2] = (int) (r / 4 + g / 2 + b / 4 +1/2); //Y

            iAux->cor[y][x].val[1] = (int) (-r / 4 + g / 2 - b / 4 + halfH + 1 / 2); //Cg
            iAux->cor[y][x].val[0] = (int) (r / 2 - b / 2 + halfH + 1 / 2); //Co
        }
    }

    // Transforma de volta para RGB

    for (y = 0; y < cI->ny; y++) {
        for (x = 0; x < cI->nx; x++) {
            Yi = iAux->cor[y][x].val[2];
            Cg = iAux->cor[y][x].val[1];
            Co = iAux->cor[y][x].val[0];

            r = Yi - Cg + Co;
            g = Yi + Cg - halfH;
            b = Yi - Cg - Co + Imax;

            if (r < 0) r = 0;
            if (g < 0) g = 0;
            if (b < 0) b = 0;

            if (r > Imax) r = Imax;
            if (g > Imax) g = Imax;
            if (b > Imax) b = Imax;

            iFinal->cor[y][x].val[2] = r;//R
            iFinal->cor[y][x].val[1] = g;//G
            iFinal->cor[y][x].val[0] = b;//B
        }
    }

    WriteColorImage( iFinal, "juncao.pgm");
}

void colorlabel( char *filename, char *filenameIN ){
    int i;
    int Imax=0, Imax2=0;
    int x, y;
    float V;
    float H = 255;
    float aux, ratio;
    int aux1, aux2, aux3, limiar=1000;

    GrayImage *iLabel = NULL;
    GrayImage *iGray = NULL;
    ColorImage *cI = NULL;


    i = 100;

    getslice( filename, "z", i, "x", "y", "00",  "brain.pgm" );
    getslice( filenameIN, "z", i, "x", "y", "00",  "brain_label.pgm" );

    iLabel = ReadGrayImage("brain_label.pgm");
    iGray = ReadGrayImage("brain.pgm");

    for( y=0; y < iGray->ny; y++ ){
        for( x=0; x < iGray->nx; x++ ){
            if (iGray->val[y][x]>limiar)
                iGray->val[y][x] = limiar;
        }
    }

    ratio = (float) 255/limiar;
    for( y=0; y < iGray->ny; y++ ){
        for( x=0; x < iGray->nx; x++ ){
            iGray->val[y][x] = iGray->val[y][x]*ratio;
        }
    }

    for( y=0; y < iLabel->ny; y++ ){
        for( x=0; x < iLabel->nx; x++ ){
            if( iLabel->val[y][x] > Imax )
                Imax = iLabel->val[y][x];
        }
    }

    if( Imax > 0 ){
        cI = CreateColorImage( iLabel->nx, iLabel->ny );

        for( y = 0; y < iLabel->ny; y++ ){
            for( x = 0; x < iLabel->nx; x++ ){
                if( iLabel->val[y][x] != 0  ){
                    aux = iLabel->val[y][x];
                    V = aux/Imax;
                    V = 4*V + 1;

                    aux1 = (int) ( H * ( Maximo( 0, ( 3 - abs( V - 4 ) - abs( V - 5 ) ) / 2 ) ) ); // Red
                    aux2 = (int) ( H * ( Maximo( 0, ( 4 - abs( V - 2 ) - abs( V - 4 ) ) / 2 ) ) ); // Green
                    aux3 = (int) ( H * ( Maximo( 0, ( 3 - abs( V - 1 ) - abs( V - 2 ) ) / 2 ) ) ); // Blue


                    cI->cor[y][x].val[2] = aux1;
                    cI->cor[y][x].val[1] = aux2;
                    cI->cor[y][x].val[0] = aux3;
                }
                else{
                    cI->cor[y][x].val[2] = 0;
                    cI->cor[y][x].val[1] = 0;
                    cI->cor[y][x].val[0] = 0;
                }
            }
        }
        WriteColorImage(cI, "cor_labels.pgm");
        composicaoImagem( iGray, cI );
    }

}

int main (){

    colorlabel("/home/rodolfo/MO815/Projeto/images/brain.scn", "/home/rodolfo/MO815/Projeto/images/brain_label.scn");

    return 0;
}
