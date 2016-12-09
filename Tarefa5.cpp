#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "mc920.h"
#include "iostream"

struct voxel{ // Voxels da imagem criada
    Matrix *coord;
    int Int;
};

MedicalImage *MedImage;
GrayImage* GetSlice(Point ponto, Vector Vet);
int cont = 0;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_OpenImage_clicked()
{
    int Imax = 0;
    MedImage = ReadMedicalImage("/home/rodolfo/MO815/Projeto/images/brain.scn");
    ui->Status->setText("Starting...");
    for (int z=0; z < MedImage->nz; z++)
        for (int y=0; y < MedImage->ny; y++)
            for (int x=0; x < MedImage->nx; x++)
                if (MedImage->val[z][y][x] >  Imax)
                    Imax = MedImage->val[z][y][x];
    MedImage->Imax=Imax;
    Point p1;
// Define os pontos iniciais
    p1.x = MedImage->nx/2;
    p1.y = MedImage->ny/2;
    p1.z = MedImage->nz/2;
    Vector Vet;
    Vet.x = -10;
    Vet.y = 0;
    Vet.z = 0;
    int aux1 = -10;
    int aux2 = 0;

// Chama várias vezes a função mudando o vetor de visualização
    for (int roda=0; roda<10; roda++) {
        GrayImage *SliceCinza = GetSlice(p1, Vet);
        Vet.x = aux1+1;
        Vet.y = aux2-1;
        aux1++;
        aux2--;
        float raiz = sqrt(Vet.x*Vet.x + Vet.y*Vet.y + Vet.z*Vet.z);
        Vet.x = Vet.x/raiz;
        Vet.y = Vet.y/raiz;
        Vet.z = Vet.z/raiz;
        printf("x:%f  y:%f z:%f\n", Vet.x, Vet.y, Vet.z);
    }
    ui->Status->setText("Done!");
}

GrayImage* GetSlice(Point ponto, Vector Vet){
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
    Matrix *rotz = RotationMatrix(AXIS_Z, (-1)*az*180/PI);
    printf("---ax:%f ay:%f az:%f\n", ax ,ay);

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
            matrix_aux = MultMatrices(rotx, matrix_aux);
            matrix_aux = MultMatrices(roty, matrix_aux);
            matrix_aux = MultMatrices(rotz, matrix_aux);
            matrix_aux = MultMatrices(transl_p1, matrix_aux);

            if ((matrix_aux->val[0] < 0) || (matrix_aux->val[1] < 0) || (matrix_aux->val[2] < 0)
                    || (matrix_aux->val[0] >= MedImage->nx) || (matrix_aux->val[1] >= MedImage->ny)
                    || (matrix_aux->val[2] >= MedImage->nz))
                SliceCinza->val[j][i] = 0;
            else {
                ponto_aux.x = matrix_aux->val[0];
                ponto_aux.y = matrix_aux->val[1];
                ponto_aux.z = matrix_aux->val[2];
                if (MedImage->Imax > 255)
                    SliceCinza->val[j][i] = (int) (ImageValueAtPoint(MedImage, ponto_aux)*255/MedImage->Imax);
                else
                    SliceCinza->val[j][i] = ImageValueAtPoint(MedImage, ponto_aux);
                DestroyMatrix(&matrix_aux);
            }
        }
    }

// Salva o arquivo
    char filename[20] = "Arquivo ";
    sprintf(filename, "%s 1%d",filename, cont);
    WriteGrayImage(SliceCinza, filename);
    cont++;

    return SliceCinza;

}


