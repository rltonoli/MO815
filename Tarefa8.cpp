#include "mainwindow.h"
#include <QApplication>
#include "mc920.h"
#include <string>
#include <time.h>

MedicalImage *MedImage;
GrayImage* Render(MedicalImage *MedImage, int AnguloX, int AnguloY, int AnguloZ, Vector VetorNormal, float ka, float kd, float ks, int ns, int La, int Lf, int raio);
int sign( int x );

struct voxel{ // Voxel do Plano de Visualização
    Matrix *coord; // Coordenadas no plano
    Matrix *pi_coord; // Coordenadas transformadas por phi_inv
    int Int; // Valor da intensidade
    float Dist;
    Vector G; // Vetor Normal
    float PhongTht;
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
    Ttx = 90;
    Tty = 0;
    Ttz = 0;
    Vet.x = 0;
    Vet.y = 0;
    Vet.z = 1;

    // Lê a imagem médica e normaliza se precisar
    MedImage = ReadMedicalImage("/home/rodolfo/MO815/Projeto/images/skullLabel60.scn");
    // Se  quiser tirar várias imagens ao mesmo tempo tem esse loop
    int final=1;
    for (int i=0; i<final; i++){

        clock_t begin1 = clock();

        // Chama o Render
        GrayImage *ImgRend = Render(MedImage, Ttx, Tty, Ttz,Vet, 0.1, 0.4, 0.5, 5, 55, 200, 20);

        // Monta o nome do arquivo para salvar
        char filename[10] = "Render";
        sprintf(filename, "%s%d-%d.pgm",filename, Ttx, Tty);
        WriteGrayImage(ImgRend, filename);
        printf("Arquivo %s salvo\n", filename);

        // Incrementa o ângulo para a próxima vez
        Tty += 10;
        DestroyGrayImage(&ImgRend);

        // Verifica quanto tempo durou a execução
        clock_t end1 = clock();
        double time_spent1 = (double) (end1-begin1)/ CLOCKS_PER_SEC;
        printf("Tempo: %f\n", time_spent1);
    }

    printf("\nTerminou\n");
    return 0;
}

GrayImage* Render(MedicalImage *MedImage, int AnguloX, int AnguloY, int AnguloZ, Vector VetorNormal, float ka, float kd, float ks, int ns, int La, int Lf, int raio){
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
    int diagonal = (sqrt(((MedImage->nx-1)*(MedImage->nx-1))+((MedImage->ny-1)*(MedImage->ny-1))+((MedImage->nz-1)*(MedImage->nz-1)) - 1));
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
            PlanoVis[i+j*diagonal].Dist = 0;
            PlanoVis[i+j*diagonal].G.x = 0;
            PlanoVis[i+j*diagonal].G.y = 0;
            PlanoVis[i+j*diagonal].G.z = 0;
            PlanoVis[i+j*diagonal].PhongTht = 2*3.145159;
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
    int Neighbor[9], count2max=0, xVisitante, yVisitante, zVisitante, Visitante, IntThreshold;
    int xaux[9], yaux[9], zaux[9];
    float termo1, termo2, termo3, lambda[6], p1[3], p2[3], lmax, lmin, MaxDist = 0,  distx, disty, distz, modvec;
    AdjRel *Adj = Spheric(raio);
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
            matrix_aux->val[3] = 0; // zerando a quarta dimensão para não zoar o produto interno
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
    //      normalmente, esse "ruído" deve ser de alguma aproximação do lambda. No limite de cima também acontecia isso,
    //      o limite é 200 mas um dos pontos ficava em 200.00001 e caía pra fora, o que está errado.
    //      Mas o arredondamento que estou fazendo está na primeira casa decimal, está muito forte, talvez um arredondamento
    //      lá pela terceira casa ficaria melhor.
                achaP->val[0] = round((PlanoVis[i+j*diagonal].pi_coord->val[0] + lambda[cont]*pi_VetorNormVis->val[0])*1000)/1000;
                achaP->val[1] = round((PlanoVis[i+j*diagonal].pi_coord->val[1] + lambda[cont]*pi_VetorNormVis->val[1])*1000)/1000;
                achaP->val[2] = round((PlanoVis[i+j*diagonal].pi_coord->val[2] + lambda[cont]*pi_VetorNormVis->val[2])*1000)/1000;
                achaP->val[3] = 1;
                if ((achaP->val[0] < MedImage->nx) && (achaP->val[1] < MedImage->ny) && (achaP->val[2] < MedImage->nz) &&
                        (achaP->val[0] >= 0) && (achaP->val[1] >= 0) && (achaP->val[2] >= 0)) {
                   //printf("%f %f %f\n", achaP->val[0], achaP->val[1], achaP->val[2]);
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
//                printf("P1: %f %f %f\n", p1[0],p1[1],p1[2]);
//                printf("PN: %f %f %f\n", p2[0],p2[1],p2[2]);
                if (n>1){
                    Point P;
                    int Ip=0;
                    IntThreshold=0;
                    P.x = p1[0];
                    P.y = p1[1];
                    P.z = p1[2];
                    for (int k = 0; k < n; k++){
                        Neighbor[0]=0;Neighbor[1]=0;Neighbor[2]=0;Neighbor[3]=0;Neighbor[4]=0;Neighbor[5]=0;Neighbor[6]=0;Neighbor[7]=0;Neighbor[8]=0;
                        //Verifica quem são os vizinhos
                        xaux[1] = (int) ceil(P.x); yaux[1] = (int) ceil(P.y); zaux[1] = (int) ceil(P.z);
                        if ((xaux[1]<MedImage->nx) && (yaux[1]<MedImage->ny) && (zaux[1]<MedImage->nz))
                            Neighbor[1] = MedImage->val[zaux[1]][yaux[1]][xaux[1]];
                        xaux[2] = (int) ceil(P.x); yaux[2] = (int) ceil(P.y); zaux[2] = (int) floor(P.z);
                        if ((xaux[2]<MedImage->nx) && (yaux[2]<MedImage->ny) && (zaux[2]<MedImage->nz))
                            Neighbor[2] = MedImage->val[zaux[2]][yaux[2]][xaux[2]];
                        xaux[3] = (int) ceil(P.x); yaux[3] = (int) floor(P.y); zaux[3] = (int) floor(P.z);
                        if ((xaux[3]<MedImage->nx) && (yaux[3]<MedImage->ny) && (zaux[3]<MedImage->nz))
                            Neighbor[3] = MedImage->val[zaux[3]][yaux[3]][xaux[3]];
                        xaux[4] = (int) ceil(P.x); yaux[4] = (int) floor(P.y); zaux[4] = (int) ceil(P.z);
                        if ((xaux[4]<MedImage->nx) && (yaux[4]<MedImage->ny) && (zaux[4]<MedImage->nz))
                            Neighbor[4] = MedImage->val[zaux[4]][yaux[4]][xaux[4]];
                        xaux[5] = (int) floor(P.x); yaux[5] = (int) ceil(P.y); zaux[5] = (int) ceil(P.z);
                        if ((xaux[5]<MedImage->nx) && (yaux[5]<MedImage->ny) && (zaux[5]<MedImage->nz))
                            Neighbor[5] = MedImage->val[zaux[5]][yaux[5]][xaux[5]];
                        xaux[6] = (int) floor(P.x); yaux[6] = (int) floor(P.y); zaux[6] = (int) ceil(P.z);
                        if ((xaux[6]<MedImage->nx) && (yaux[6]<MedImage->ny) && (zaux[6]<MedImage->nz))
                            Neighbor[6] = MedImage->val[zaux[6]][yaux[6]][xaux[6]];
                        xaux[7] = (int) floor(P.x); yaux[6] = (int) ceil(P.y); zaux[7] = (int) floor(P.z);
                        if ((xaux[7]<MedImage->nx) && (yaux[7]<MedImage->ny) && (zaux[7]<MedImage->nz))
                            Neighbor[7] = MedImage->val[zaux[7]][yaux[7]][xaux[7]];
                        xaux[8] = (int) floor(P.x); yaux[8] = (int) floor(P.y); zaux[8] = (int) floor(P.z);
                        if ((xaux[8]<MedImage->nx) && (yaux[8]<MedImage->ny) && (zaux[8]<MedImage->nz))
                            Neighbor[8] = MedImage->val[zaux[8]][yaux[8]][xaux[8]];
                        Ip = 0;
			// Pega o valor máximo dos vizinhos
                        for (int count2=0; count2 <9; count2++){
                            if (Neighbor[count2] > Ip) {
                                Ip = Neighbor[count2];
                                count2max = count2;
                            }
                        }
                        // Verifica se é um pixel de suerfície (daria pra verificar se é de label aqui na tarefa de transparencia)
                        if (Ip > IntThreshold)  {
                            // Calcula Depth Shading
                            PlanoVis[i+j*diagonal].Int = ka*La;
                            distx = PlanoVis[i+j*diagonal].pi_coord->val[0]-xaux[count2max];
                            disty = PlanoVis[i+j*diagonal].pi_coord->val[1]-yaux[count2max];
                            distz = PlanoVis[i+j*diagonal].pi_coord->val[2]-zaux[count2max];
                            PlanoVis[i+j*diagonal].Dist = sqrt(xaux[count2max]*xaux[count2max] + yaux[count2max]*yaux[count2max] + zaux[count2max]*zaux[count2max]);
                            if (PlanoVis[i+j*diagonal].Dist > MaxDist) MaxDist = PlanoVis[i+j*diagonal].Dist;

                            // Estima o vetor normal para o calculo da iluminação de phong
                            for (int iadj=0;iadj<Adj->n;iadj++){
				// Visita os pixels da adjacencia esferica
                                xVisitante = xaux[count2max]+Adj->adj[iadj].dx;
                                yVisitante = yaux[count2max]+Adj->adj[iadj].dy;
                                zVisitante = zaux[count2max]+Adj->adj[iadj].dz;
                                if ((xVisitante < MedImage->nx) && (yVisitante < MedImage->ny) && (zVisitante < MedImage->nz)
                                        &&  (xVisitante >= 0) &&  (yVisitante >= 0) &&  (zVisitante >= 0)) {
                                    Visitante = MedImage->val[zVisitante][yVisitante][xVisitante];
                                    distx = xaux[count2max] - xVisitante;
                                    disty = yaux[count2max] - yVisitante;
                                    distz = zaux[count2max] - zVisitante;
                                    if (Visitante > IntThreshold) {// && faz parte do objeto na transparencia
                                        modvec = sqrt(distx*distx + disty*disty + distz*distz);
                                        if (modvec!=0){
                                            PlanoVis[i+j*diagonal].G.x += Visitante*distx/modvec;
                                            PlanoVis[i+j*diagonal].G.y += Visitante*disty/modvec;
                                            PlanoVis[i+j*diagonal].G.z += Visitante*distz/modvec;
                                        }
                                    }
                                }
                            }
                            if ((PlanoVis[i+j*diagonal].G.x!=0) || (PlanoVis[i+j*diagonal].G.y!=0) || (PlanoVis[i+j*diagonal].G.z!=0)){
				// Calcula o modulo do vetor
                                modvec = sqrt(pow(pi_VetorNormVis->val[0],2)+pow(pi_VetorNormVis->val[1],2)+pow(pi_VetorNormVis->val[2],2))*sqrt(pow(PlanoVis[i+j*diagonal].G.x,2)+pow(PlanoVis[i+j*diagonal].G.y,2)+pow(PlanoVis[i+j*diagonal].G.z,2));
				// Salva o valor do angulo do vetor de visualização até o vetor normal da superfície
                                PlanoVis[i+j*diagonal].PhongTht = acos((pi_VetorNormVis->val[0]*PlanoVis[i+j*diagonal].G.x+pi_VetorNormVis->val[1]*PlanoVis[i+j*diagonal].G.y+pi_VetorNormVis->val[2]*PlanoVis[i+j*diagonal].G.z)/modvec);
                            }
                            k = n; // Para de andar na reta do DDA já que só interessa a superficie (na transparencia continua)
                        }
                        // Da o próximo passo na reta
                        P.x += dx;
                        P.y += dy;
                        P.z += dz;
                    }
                }
            }

        }
    }
// Joga para uma imagem cinza calculando o modelo de phong
    GrayImage *GrayImg = CreateGrayImage(diagonal,diagonal);
    for (int j=0;j<diagonal;j++){
        for (int i=0;i<diagonal;i++){
            PlanoVis[i+j*diagonal].Dist = PlanoVis[i+j*diagonal].Dist/MaxDist;
            GrayImg->val[j][i] = PlanoVis[i+j*diagonal].Int;
            if ((PlanoVis[i+j*diagonal].PhongTht < 3.14159) && (PlanoVis[i+j*diagonal].PhongTht > 3.14159/2)) GrayImg->val[j][i] += (int) abs(Lf*kd*cos(PlanoVis[i+j*diagonal].PhongTht)*(1-PlanoVis[i+j*diagonal].Dist));
            if ((PlanoVis[i+j*diagonal].PhongTht > 3*3.14159/4) && (PlanoVis[i+j*diagonal].PhongTht <3.14159)) GrayImg->val[j][i] += (int) abs(Lf*ks*pow(cos(2*PlanoVis[i+j*diagonal].PhongTht),ns)*(1-PlanoVis[i+j*diagonal].Dist));

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
