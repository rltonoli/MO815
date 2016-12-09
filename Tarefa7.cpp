#include "mainwindow.h"
#include <QApplication>
#include "mc920.h"
#include <string>
#include <time.h>

MedicalImage *MedImage;
GrayImage* GetMIP(MedicalImage *MedImage, int AnguloX, int AnguloY, int AnguloZ,Vector VetorNormal, int modo);
int sign( int x );
MedicalImage* Normaliza(MedicalImage *MedImage);

struct voxel{ // Voxel do Plano de Visualização
    Matrix *coord; // Coordenadas no plano
    Matrix *pi_coord; // Coordenadas transformadas por phi_inv
    int Int; // Valor da intensidade
};

struct Faces{ // Struct contendo informações de cada face (6 no total)
    Matrix *VetNorm; // Vetor Normal da Face
    Matrix *Centro; // Centro da Face
};

int main(int argc, char *argv[])
{
    // Define os angulos de rotação e seta o vetor de visualização (no caso do MIP ele é fixo)
    int Ttx, Tty, Ttz;
    Vector Vet;
    int Modo;
    Ttx = 80;
    Tty = -10;
    Ttz = 0;
    Vet.x = 0;
    Vet.y = 0;
    Vet.z = 1;
    Modo = 1;

    // Lê a imagem médica
    MedImage = ReadMedicalImage("/home/rodolfo/MO815/Projeto/images/brain.scn");

    // Se  quiser tirar várias imagens ao mesmo tempo tem esse loop
    int final=3;
    for (int i=0; i<final; i++){

        clock_t begin1 = clock();

        // Chama o MIP
        GrayImage *MIP = GetMIP(MedImage, Ttx, Tty, Ttz,Vet, Modo);

        // Monta o nome do arquivo para salvar
        if (Modo == 0) {
            char filename[10] = "MIP";
            sprintf(filename, "%s%d.pgm",filename, Tty);
            WriteGrayImage(MIP, filename);
            printf("Arquivo %s salvo\n", filename);
        }

        if (Modo == 1) {
            char filename[13] = "MeanIP";
            sprintf(filename, "%s%d.pgm",filename, Tty);
            WriteGrayImage(MIP, filename);
            printf("Arquivo %s salvo\n", filename);
        }

        // Incrementa o ângulo para a próxima vez
        Tty += 10;
        DestroyGrayImage(&MIP);

        // Verifica quanto tempo durou a execução
        clock_t end1 = clock();
        double time_spent1 = (double) (end1-begin1)/ CLOCKS_PER_SEC;
        printf("Tempo: %f\n", time_spent1);
    }

    printf("\nTerminou\n");
    return 0;
}

GrayImage* GetMIP(MedicalImage *MedImage, int AnguloX, int AnguloY, int AnguloZ, Vector VetorNormal, int modo){
// MODO 0 = Projeção Maxima. MODO = 1 Projeção Media
// --- Criação das matrizes de translação
    Matrix *transl_pc = CreateMatrix(4,4);
    Matrix *transl_qc = CreateMatrix(4,4);

    // Identidade de translação de pc e qc
    for (int i=0;i<4;i++) {
        for (int j=0; j<4;j++) {
            if (i==j) {
                transl_pc->val[GetMatrixIndex(transl_pc,i,j)] = 1;
                transl_qc->val[GetMatrixIndex(transl_qc,i,j)] = 1;
            }
            else if ((i!=3) && (i!=j)) {
                transl_pc->val[GetMatrixIndex(transl_pc,i,j)] = 0;
                transl_qc->val[GetMatrixIndex(transl_qc,i,j)] = 0;
            }
        }
    }

    // Valores de translação de pc e qc
    // Estou tirando um pois o último pixel pq o MedImage->nx,ny e nz são tamanhos e não posições.
    // Posições contam o zero e vão de 0 até 99 com o tamanho 100, por exemplo
    int diagonal = (sqrt(((MedImage->nx-1)*(MedImage->nx-1))+((MedImage->ny-1)*(MedImage->ny-1))+((MedImage->nz-1)*(MedImage->nz-1)) + 1));
    transl_pc->val[GetMatrixIndex(transl_pc,3,0)] = (int) (MedImage->nx-1)/2;
    transl_pc->val[GetMatrixIndex(transl_pc,3,1)] = (int) (MedImage->ny-1)/2;
    transl_pc->val[GetMatrixIndex(transl_pc,3,2)] = (int) (MedImage->nz-1)/2;
    transl_qc->val[GetMatrixIndex(transl_qc,3,0)] = -(int) (diagonal/2);
    transl_qc->val[GetMatrixIndex(transl_qc,3,1)] = -(int) (diagonal/2);
    transl_qc->val[GetMatrixIndex(transl_qc,3,2)] = -(int) (diagonal/2);

    // Criação das matrizes de rotação
    Matrix *rotx = RotationMatrix(AXIS_X,(-1)*AnguloX);
    Matrix *roty = RotationMatrix(AXIS_Y,(-1)*AnguloY);
    Matrix *rotz = RotationMatrix(AXIS_Z,(-1)*AnguloZ);

// --- Criação do plano de visualização (q). Cada elemento x do PlanoVis[x] é um pixel do plano de visualização.
    // Cada elemento possui uma matriz correspondente a suas coordenadas normais e outra com elas transformadas por i_inv.
    struct voxel PlanoVis[diagonal*diagonal];
    for (int j=0;j<diagonal;j++){
        for (int i=0;i<diagonal;i++){
            PlanoVis[i+j*diagonal].coord = CreateMatrix(1,4);
            PlanoVis[i+j*diagonal].coord->val[GetMatrixIndex(PlanoVis[i+j*diagonal].coord,0,0)] = i;
            PlanoVis[i+j*diagonal].coord->val[GetMatrixIndex(PlanoVis[i+j*diagonal].coord,0,1)] = j;
            PlanoVis[i+j*diagonal].coord->val[GetMatrixIndex(PlanoVis[i+j*diagonal].coord,0,2)] = -(int) (diagonal/2);
            PlanoVis[i+j*diagonal].coord->val[GetMatrixIndex(PlanoVis[i+j*diagonal].coord,0,3)] = 1;
            PlanoVis[i+j*diagonal].Int = 0;
        }
    }

// --- Criaçao dos valores das faces de um cubo
    // Serão criadas como matrizes 1x4 pq PlanoVis e VetNorm também são 1x4, mas a última dimensão
    // será zerada para não afetar em nada
    struct Faces Face[6];
    for (int cont = 0; cont < 6; cont++) {
        Face[cont].VetNorm = CreateMatrix(1,4);
        Face[cont].Centro = CreateMatrix(1,4);
        Face[cont].VetNorm->val[3] = 0;
        Face[cont].Centro->val[3] = 0;
        switch (cont) {
        case 0: // Frente
            Face[0].VetNorm->val[0] = 0;
            Face[0].VetNorm->val[1] = 0;
            Face[0].VetNorm->val[2] = -1;
            Face[0].Centro->val[0] = (MedImage->nx-1)/2;
            Face[0].Centro->val[1] = (MedImage->ny-1)/2;
            Face[0].Centro->val[2] = 0;
           break;
        case 1: // Costas
            Face[1].VetNorm->val[0] = 0;
            Face[1].VetNorm->val[1] = 0;
            Face[1].VetNorm->val[2] = 1;
            Face[1].Centro->val[0] = (MedImage->nx-1)/2;
            Face[1].Centro->val[1] = (MedImage->ny-1)/2;
            Face[1].Centro->val[2] = MedImage->nz-1;
            break;
        case 2: // Direita
            Face[2].VetNorm->val[0] = 1;
            Face[2].VetNorm->val[1] = 0;
            Face[2].VetNorm->val[2] = 0;
            Face[2].Centro->val[0] = MedImage->nx-1;
            Face[2].Centro->val[1] = (MedImage->ny-1)/2;
            Face[2].Centro->val[2] = (MedImage->nz-1)/2;
            break;
        case 3: // Esquerda
            Face[3].VetNorm->val[0] = -1;
            Face[3].VetNorm->val[1] = 0;
            Face[3].VetNorm->val[2] = 0;
            Face[3].Centro->val[0] = 0;
            Face[3].Centro->val[1] = (MedImage->ny-1)/2;
            Face[3].Centro->val[2] = (MedImage->nz-1)/2;
            break;
        case 4: // Cima
            Face[4].VetNorm->val[0] = 0;
            Face[4].VetNorm->val[1] = -1;
            Face[4].VetNorm->val[2] = 0;
            Face[4].Centro->val[0] = (MedImage->nx-1)/2;
            Face[4].Centro->val[1] = 0;
            Face[4].Centro->val[2] = (MedImage->nz-1)/2;
            break;
        case 5: // Baixo
            Face[5].VetNorm->val[0] = 0;
            Face[5].VetNorm->val[1] = 1;
            Face[5].VetNorm->val[2] = 0;
            Face[5].Centro->val[0] = (MedImage->nx-1)/2;
            Face[5].Centro->val[1] = MedImage->ny-1;
            Face[5].Centro->val[2] = (MedImage->nz-1)/2;
            break;
        default:
            break;
        }
    }

// --- Calcula phi_inv do vetor normal de visualização
    Matrix *VetorNormVis = CreateMatrix(1,4);
    Matrix *pi_VetorNormVis = CreateMatrix(1,4);
    VetorNormVis->val[0] = VetorNormal.x;
    VetorNormVis->val[1] = VetorNormal.y;
    VetorNormVis->val[2] = VetorNormal.z;
    VetorNormVis->val[3] = 1;

// O vetor de visualização não precisa ser transladado pois ele só dá a direção
    //pi_VetorNormVis = MultMatrices(transl_qc, VetorNormVis);
    pi_VetorNormVis = MultMatrices(rotz, VetorNormVis);
    pi_VetorNormVis = MultMatrices(roty, pi_VetorNormVis);
    pi_VetorNormVis = MultMatrices(rotx, pi_VetorNormVis);
    //pi_VetorNormVis = MultMatrices(transl_pc, pi_VetorNormVis);
    pi_VetorNormVis->val[3] = 0;
    DestroyMatrix(&VetorNormVis);

// --- Aqui é onde a maior parte das coisas acontecem
// --- 1: Cada loop é um pixel do Plano de Visualização (ponto Q)
// --- 2: Paro Q' = phi_inv(Q))
// --- 3: Para cada face é calculado o ponto de intersecção e verificado se está dentro da cena
// --- 4: Caso positivo calcula-se P = Q' + lambda * n' para descobrir P1 e PN
// --- 5: Aplica DDA em três dimensões guardando o valor máximo da reta
    Matrix *matrix_aux = CreateMatrix(1,4);
    Matrix *matrix_aux2 = CreateMatrix(1,1);
    Matrix *achaP = CreateMatrix(1,4);
    // Para realizar a rotação os vetores precisam ter a quarta dimensão com valor 1, mas isso da erro no calculo
    // do produto interno, então essa dimensão será zerada
    int cont_debug=0;
    float termo1, termo2, termo3, lambda[6], p1[3], p2[3], lmax, lmin;
// 1:
    for (int j=0;j<diagonal;j++){
        for (int i=0;i<diagonal;i++){
            lmax = 0;
            lmin = 2*diagonal;
// 2:
            matrix_aux = MultMatrices(transl_qc, PlanoVis[i+j*diagonal].coord);
            matrix_aux = MultMatrices(rotz, matrix_aux);
            matrix_aux = MultMatrices(roty, matrix_aux);
            matrix_aux = MultMatrices(rotx, matrix_aux);
            matrix_aux = MultMatrices(transl_pc, matrix_aux);
            matrix_aux->val[3] = 0; // zerando a quarta dimensão para não dar erro no produto interno
            PlanoVis[i+j*diagonal].pi_coord = CreateMatrix(1,4);
            PlanoVis[i+j*diagonal].pi_coord->val[GetMatrixIndex(PlanoVis[i+j*diagonal].coord,0,0)] = matrix_aux->val[0];
            PlanoVis[i+j*diagonal].pi_coord->val[GetMatrixIndex(PlanoVis[i+j*diagonal].coord,0,1)] = matrix_aux->val[1];
            PlanoVis[i+j*diagonal].pi_coord->val[GetMatrixIndex(PlanoVis[i+j*diagonal].coord,0,2)] = matrix_aux->val[2];
            PlanoVis[i+j*diagonal].pi_coord->val[GetMatrixIndex(PlanoVis[i+j*diagonal].coord,0,3)] = matrix_aux->val[3];
            DestroyMatrix(&matrix_aux);
// 3:
            for (int cont=0; cont<6; cont++){ // ( nj*cj/mod - nj*phi_inv(q)/mod ) / nj*phi_inv(n)/mod
                matrix_aux2 = MultMatrices(TransposeMatrix(Face[cont].VetNorm), (Face[cont].Centro));
                termo1 = matrix_aux2->val[0];
                matrix_aux2 = MultMatrices(TransposeMatrix(Face[cont].VetNorm), PlanoVis[i+j*diagonal].pi_coord);
                termo2 = matrix_aux2->val[0];
                matrix_aux2 = MultMatrices(TransposeMatrix(Face[cont].VetNorm), pi_VetorNormVis);
                termo3 = matrix_aux2->val[0];
                lambda[cont] = (termo1-termo2)/termo3;

    // !!!  Aqui eu faço o "round" pq vários pontos davam QUASE no limite da cena mas eram deixados de fora, por exemplo
    //      -0.000003 ficava de fora do calculo de lambda, mas na verdade ele era para ser 0 redondo e calcular o lambda
    //      normalmente, esse erro deve ser de alguma aproximação na hora do calculo da rotação. No limite de cima também acontecia isso,
    //      por exemplo o limite é 200 mas um dos pontos ficava em 200.00001 e caia pra fora, o que não era para acontecer.
                achaP->val[0] = round((PlanoVis[i+j*diagonal].pi_coord->val[0] + lambda[cont]*pi_VetorNormVis->val[0])*1000)/1000;
                achaP->val[1] = round((PlanoVis[i+j*diagonal].pi_coord->val[1] + lambda[cont]*pi_VetorNormVis->val[1])*1000)/1000;
                achaP->val[2] = round((PlanoVis[i+j*diagonal].pi_coord->val[2] + lambda[cont]*pi_VetorNormVis->val[2])*1000)/1000;
                achaP->val[3] = 1;
                if ((achaP->val[0] < MedImage->nx) && (achaP->val[1] < MedImage->ny) && (achaP->val[2] < MedImage->nz) &&
                        (achaP->val[0] >= 0) && (achaP->val[1] >= 0) && (achaP->val[2] >= 0)) {
// 4:
                    cont_debug++; // Contador usado para ver quantos pontos válidos foram encontrados
                    // O maior lambda válido é o mais distante, o ponto PN
                    if (lambda[cont] > lmax) {
                        lmax = lambda[cont];
                        p2[0] = achaP->val[0];
                        p2[1] = achaP->val[1];
                        p2[2] = achaP->val[2];
                    }
                    // O menor lambda válido é o que intercepta o plano mais próximo, o ponto P1
                    if (lambda[cont] < lmin) {
                        lmin = lambda[cont];
                        p1[0] = achaP->val[0];
                        p1[1] = achaP->val[1];
                        p1[2] = achaP->val[2];
                    }
                }
            }
// 5:
            float Dx, Dy, Dz, dx, dy, dz;
            int n;
            if ((p1[0] == p2[0]) && (p1[1] == p2[1]) && (p1[2] == p2[2])) {
                // Se cair aqui em algum momento, algo está errado
                n = 1;
            }
            else {
                Dx = p2[0] - p1[0];
                Dy = p2[1] - p1[1];
                Dz = p2[2] - p1[2];
                if( ( abs(Dx) >= abs(Dy) ) && ( abs(Dx) >= abs(Dz) ) ){
                    n = abs(Dx)+1;
                    dx = sign(Dx);
                    dy = dx * (Dy/Dx);
                    dz = dx * (Dz/Dx);
                }
                else{
                    if( ( abs(Dy) >= abs(Dx) ) && ( abs(Dy) >= abs(Dz) ) ){
                        n = abs(Dy)+1;
                        dy = sign(Dy);
                        dx = dy * (Dx/Dy);
                        dz = dy * (Dz/Dy);
                    }
                    else{
                        n = abs(Dz)+1;
                        dz = sign(Dz);
                        dx = dz * (Dx/Dz);
                        dy = dz * (Dy/Dz);
                    }

                }
                if (n>1){
                    Point P;
                    float IpFix = 0, Ip=0;
                    int NPixelValido=0;
                    P.x = p1[0];
                    P.y = p1[1];
                    P.z = p1[2];
		    // Caminha pela reta
                    for (int k = 0; k < n; k++){
			// Verifica se pertence a cena
                        if ((P.x < MedImage->nx) && (P.y < MedImage->ny) && (P.z < MedImage->nz)){
                            Ip = ImageValueAtPoint(MedImage, P);
                            if (modo == 0){ // Projeção de intensidade máxima
                                if (Ip > IpFix) IpFix = Ip;
                            }
                            else { // Projeção de intensidade média
                                IpFix += Ip;
                                NPixelValido++;
                            }
                        }
                        P.x += dx;
                        P.y += dy;
                        P.z += dz;
                    }
                    // Joga o valor máximo ou médio no Plano de Visualização
                    if ((modo == 1) && (NPixelValido > 0)) IpFix = IpFix/NPixelValido;
                    PlanoVis[i+j*diagonal].Int = (int) IpFix;
                }
            }

        }
    }

// --- Apenas joga o Plano de Visualização em uma imagem cinza para salvar
    GrayImage *GrayImg = CreateGrayImage(diagonal,diagonal);
    for (int j=0;j<diagonal;j++){
        for (int i=0;i<diagonal;i++){
            GrayImg->val[j][i] = PlanoVis[i+j*diagonal].Int;
        }
    }

    DestroyMatrix(&matrix_aux);
    DestroyMatrix(&matrix_aux2);
    DestroyMatrix(&pi_VetorNormVis);
    DestroyMatrix(&transl_pc);
    DestroyMatrix(&transl_qc);
    DestroyMatrix(&rotx);
    DestroyMatrix(&roty);
    DestroyMatrix(&rotz);
    DestroyMatrix(&pi_VetorNormVis);
    return GrayImg;
}

int sign( int x ){
    if(x > 0)
        return 1;
    if(x < 0)
        return -1;
    return 0;
}


