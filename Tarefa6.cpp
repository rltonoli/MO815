#include "mainwindow.h"
#include <QApplication>
#include "mc920.h"

struct voxel{ // Voxels da imagem criada
    Matrix *coord;
    int Int;
};

MedicalImage *MedImage;
GrayImage* GetSlice(Point ponto, Vector Vet, int count);
int Reform(MedicalImage *MedImage, Point pini, Point pfim, int cortes);
int sign( int x );

int main(int argc, char *argv[])
{
// Abre a imagem médica, define os dois pontos e chama a função
    MedImage = ReadMedicalImage("/home/rodolfo/MO815/Projeto/images/brain.scn");
    Point p1;
    p1.x = 0;
    p1.y = 0;
    p1.z = 0;
    Point p2;
    p2.x = MedImage->nx;
    p2.y = MedImage->ny;
    p2.z = MedImage->nz;
    Reform(MedImage, p1, p2, 40);
}

int Reform(MedicalImage *MedImage, Point pini, Point pfim, int cortes){
    int Imax = 0;
    float raiz;
    Point P, Passo;
    GrayImage* Slices[cortes];
    // Pega o Imax de qq jeito (8 bits não pega na lib)
    for (int z=0; z < MedImage->nz; z++)
        for (int y=0; y < MedImage->ny; y++)
            for (int x=0; x < MedImage->nx; x++)
                if (MedImage->val[z][y][x] >  Imax)
                    Imax = MedImage->val[z][y][x];
    MedImage->Imax=Imax;
    Vector Vet;
// Define o vetor de visualização
    Vet.x = pfim.x - pini.x;
    Vet.y = pfim.y - pini.y;
    Vet.z = pfim.z - pini.z;
    raiz = sqrt(Vet.x*Vet.x + Vet.y*Vet.y + Vet.z*Vet.z);
// Deixa o vetor unitário
    if (raiz!=0){
        Vet.x = Vet.x/raiz;
        Vet.y = Vet.y/raiz;
        Vet.z = Vet.z/raiz;
    }
    printf("%f %f %f\n", Vet.x, Vet.y, Vet.z);
// Calcula o passo
    Passo.x = (pfim.x - pini.x)/cortes;
    Passo.y = (pfim.y - pini.y)/cortes;
    Passo.z = (pfim.z - pini.z)/cortes;
    P.x = pini.x;
    P.y = pini.y;
    P.z = pini.z;
// Chama a função da tarefa 5 quantas vezes precisar
    for (int k = 0; k < cortes; k++){
        Slices[k] = GetSlice(P, Vet, k);
        P.x += Passo.x;
        P.y += Passo.y;
        P.z += Passo.z;
    }
    return 0;
}

// Função da tarefa anterior
GrayImage* GetSlice(Point ponto, Vector Vet, int count){
    float Vx = Vet.x;
    float Vy = Vet.y;
    float Vz = Vet.z;
    float ax, ay, az;

// --- Criação das matrizes de translação
    Matrix *transl_p1 = CreateMatrix(4,4);
    Matrix *transl_qc = CreateMatrix(4,4);

    // Identidade de translação de p1 e qc
    for (int i=0;i<4;i++) {
        for (int j=0; j<4;j++) {
            if (i==j) {
                transl_p1->val[GetMatrixIndex(transl_p1,i,j)] = 1;
                transl_qc->val[GetMatrixIndex(transl_qc,i,j)] = 1;
            }
            else if ((i!=3) && (i!=j)) {
                transl_p1->val[GetMatrixIndex(transl_p1,i,j)] = 0;
                transl_qc->val[GetMatrixIndex(transl_qc,i,j)] = 0;
            }
        }
    }
    // Valores de translação de p1 e qc
    int diagonal = (sqrt((MedImage->nx*MedImage->nx)+(MedImage->ny*MedImage->ny)+(MedImage->nz*MedImage->nz)) + 1);
    transl_p1->val[GetMatrixIndex(transl_p1,3,0)] = ponto.x;
    transl_p1->val[GetMatrixIndex(transl_p1,3,1)] = ponto.y;
    transl_p1->val[GetMatrixIndex(transl_p1,3,2)] = ponto.z;
    transl_qc->val[GetMatrixIndex(transl_qc,3,0)] = -(int) (diagonal/2);
    transl_qc->val[GetMatrixIndex(transl_qc,3,1)] = -(int) (diagonal/2);
    transl_qc->val[GetMatrixIndex(transl_qc,3,2)] = 0; // Verificar se tudo bem

    // --- Criação das matrizes de rotação
    // Cálculo de alpha x(ax) e alpha y (ay)
    if (Vz>0) {
        ax = atan(Vy/Vz);
        ay = atan(Vx*cos(ax)/Vz);
    }
    else if (Vz<0) {
        ax = atan(Vy/Vz);
        if (ax!=0) ax = ax - PI;
        ay = atan((Vx*cos(ax))/Vz);
        if (ay!=0) ay = ay - PI;
    }
    else {
        if ((Vx==0) && (Vy!=0)) {
            ay = 0;
            if (Vy>0) ax = PI/2;
            else ax = -PI/2;
        }
        else if ((Vx!=0) && (Vy==0)) {
            ax = 0;
            if (Vx>0) ay = PI/2;
            else ay = -PI/2;
        }
        else { // Vx!=0 e Vy!=0
            if (Vy>0) az = -acos(Vx);
            else az = acos(Vx);
            ay = -PI/2;
        }
    }
    Matrix *rotx = RotationMatrix(AXIS_X,(-1)*ax*180/PI);
    Matrix *roty = RotationMatrix(AXIS_Y,(-1)*ay*180/PI);
    Matrix *rotz = RotationMatrix(AXIS_Z,(-1)*az*180/PI);
    // Criaçao do plano de visualizaçao
    struct voxel PlanoVis[diagonal*diagonal];
    for (int j=0;j<diagonal;j++){
        for (int i=0;i<diagonal;i++){
            PlanoVis[i+j*diagonal].coord = CreateMatrix(1,4);
            PlanoVis[i+j*diagonal].coord->val[GetMatrixIndex(PlanoVis[i+j*diagonal].coord,0,0)] = i;
            PlanoVis[i+j*diagonal].coord->val[GetMatrixIndex(PlanoVis[i+j*diagonal].coord,0,1)] = j;
            PlanoVis[i+j*diagonal].coord->val[GetMatrixIndex(PlanoVis[i+j*diagonal].coord,0,2)] = 0;
            PlanoVis[i+j*diagonal].coord->val[GetMatrixIndex(PlanoVis[i+j*diagonal].coord,0,3)] = 1;
        }
    }
    // --- Aplica phi_inv para cada ponto de qc

    Matrix *matrix_aux = CreateMatrix(1,4);
    Point ponto_aux;
    GrayImage *SliceCinza = CreateGrayImage(diagonal,diagonal);
    for (int j=0;j<diagonal;j++){
        for (int i=0;i<diagonal;i++){
            matrix_aux = MultMatrices(transl_qc, PlanoVis[i+j*diagonal].coord);


            matrix_aux = MultMatrices(rotz, matrix_aux);
            matrix_aux = MultMatrices(roty, matrix_aux);
            matrix_aux = MultMatrices(rotx, matrix_aux);

            matrix_aux = MultMatrices(transl_p1, matrix_aux);

            if ((matrix_aux->val[0] < 0) || (matrix_aux->val[1] < 0) || (matrix_aux->val[2] < 0)
                    || (matrix_aux->val[0] >= MedImage->nx) || (matrix_aux->val[1] >= MedImage->ny)
                    || (matrix_aux->val[2] >= MedImage->nz))
                SliceCinza->val[j][i] = 0;
            else {
                ponto_aux.x = matrix_aux->val[0];
                ponto_aux.y = matrix_aux->val[1];
                ponto_aux.z = matrix_aux->val[2];
                SliceCinza->val[j][i] = ImageValueAtPoint(MedImage, ponto_aux);
                DestroyMatrix(&matrix_aux);
            }
        }
    }


    char filename[20] = "Arquivo ";
    sprintf(filename, "%s %d",filename, count);
    WriteGrayImage(SliceCinza, filename);

}

int sign( int x ){
    if(x > 0)
        return 1;
    if(x < 0)
        return -1;
    return 0;
}
