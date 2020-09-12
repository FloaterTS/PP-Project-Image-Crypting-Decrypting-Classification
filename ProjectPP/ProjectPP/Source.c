#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct pixel {
	unsigned char R, G, B;
};

struct detectie {
	unsigned int x, y, latime, inaltime, sablon;
	float corr;
};

unsigned int grayscaleImage(char* nume_fisier_sursa, char* nume_fisier_destinatie, int okType)
{
	FILE *fin, *fout;
	unsigned int dim_img, latime_img, inaltime_img;
	unsigned char aux, *pRGB = (unsigned char*)malloc(3);
	int i, j;

	fin = fopen(nume_fisier_sursa, "rb");
	if (fin == NULL)
	{
		printf("Nu am gasit imaginea sursa din care citesc.\n");
		return 0;
	}
	fout = fopen(nume_fisier_destinatie, "wb+");

	fseek(fin, 2, SEEK_SET);
	fread(&dim_img, sizeof(unsigned int), 1, fin);

	fseek(fin, 18, SEEK_SET);
	fread(&latime_img, sizeof(unsigned int), 1, fin);
	fread(&inaltime_img, sizeof(unsigned int), 1, fin);

	if (okType) printf("\nTransformare grayscale..");
	fseek(fin, 0, SEEK_SET);
	unsigned char c; // copiem headerul
	for (i = 0; i < 54; i++)
	{
		fread(&c, 1, 1, fin);
		fwrite(&c, 1, 1, fout);
		fflush(fout);
	}

	int padding;
	if (latime_img % 4 != 0)
		padding = 4 - (3 * latime_img) % 4;
	else
		padding = 0;

	if (okType) printf(".");
	for (i = 0; i < inaltime_img; i++)
	{
		for (j = 0; j < latime_img; j++)
		{
			fread(pRGB, 1, 3, fin);
			aux = 0.299*pRGB[2] + 0.587*pRGB[1] + 0.114*pRGB[0];
			pRGB[0] = pRGB[1] = pRGB[2] = aux;
			fwrite(pRGB, 1, 3, fout);
			fflush(fout);
		}
		if (padding)
		{
			fread(pRGB, 1, padding, fin);
			fwrite(pRGB, 1, padding, fout);
		}
		if (okType) if (i % 128 == 0) //vizualizare progres
			printf(".");
	}
	if (okType) printf("\n");
	fclose(fin);
	fclose(fout);
	return dim_img;
}

void xorshift32_C1(unsigned int **R, unsigned int n, unsigned int seed)
{
	unsigned int r = seed;

	(*R) = (unsigned int*)malloc(n * sizeof(unsigned int));

	(*R)[0] = r;
	int k;
	for (k = 1; k < n; k++)
	{
		r = r ^ r << 13;
		r = r ^ r >> 17;
		r = r ^ r << 5;
		(*R)[k] = r;
	}
}

void pDurstenfeld(unsigned int **P, unsigned int n, unsigned int *R)
{
	unsigned int r, k, rk, aux;

	(*P) = (unsigned int*)malloc(n * sizeof(unsigned int));

	for (k = 0; k < n; k++)
		(*P)[k] = k;

	rk = 1;
	for (k = n - 1; k > 0; k--)
	{
		r = R[rk] % (k + 1);
		aux = (*P)[r];
		(*P)[r] = (*P)[k];
		(*P)[k] = aux;
		rk++;
	}
}

void inversaPerm(unsigned int **P, unsigned int n)
{
	unsigned int *invP = (unsigned int*)malloc(n * sizeof(unsigned int));
	int i;
	for (i = 0; i < n; i++)
		invP[(*P)[i]] = i;

	for (i = 0; i < n; i++)
		(*P)[i] = invP[i];
	free(invP);
}

int getSecretKeys(char *nume_fisier_secret_key, unsigned int *R0, unsigned int *SV)
{
	FILE *s;
	s = fopen(nume_fisier_secret_key, "r");
	if (s == NULL)
	{
		printf("Nu am gasit fisierul ce contine cheile secrete.\n");
		return 0;
	}
	fscanf(s, "%u", &(*R0));
	fscanf(s, "%u", &(*SV));
	fclose(s);
	return 1;
}

unsigned int loadImageLin_C2(char *nume_fisier_sursa, unsigned char **header, struct pixel **L, unsigned int *latime_img, unsigned int *inaltime_img, unsigned int *padding, int okType)
{
	unsigned int dim_img;
	int i, j;

	if (okType)
		printf("\nNume_fisier_sursa = %s \n", nume_fisier_sursa);

	//Deschidere fisier imagine
	FILE *fin;
	fin = fopen(nume_fisier_sursa, "rb");
	if (fin == NULL)
	{
		printf("Nu am gasit imaginea sursa din care citesc.\n");
		return 0;
	}

	//Dimensiune imagine
	fseek(fin, 2, SEEK_SET);
	fread(&dim_img, sizeof(unsigned int), 1, fin);
	if (okType)
		printf("Dimensiunea imaginii in octeti: %u\n", dim_img);

	//Latime si inaltime imagine
	fseek(fin, 18, SEEK_SET);
	fread(latime_img, sizeof(unsigned int), 1, fin);
	fread(inaltime_img, sizeof(unsigned int), 1, fin);
	if (okType)
		printf("Dimensiunea imaginii in pixeli (latime x inaltime): %u x %u\n", *latime_img, *inaltime_img);

	//Calculare padding
	if ((*latime_img) % 4 != 0)
		*padding = 4 - (3 * (*latime_img)) % 4;
	else
		*padding = 0;
	if (okType)
		printf("padding = %d \n", *padding);

	//Copierea header-ului
	*header = (unsigned char*)malloc(54);
	fseek(fin, 0, SEEK_SET);
	for (i = 0; i < 54; i++)
		fread(&(*header)[i], 1, 1, fin);

	//Alocare spatiu pt imaginea liniarizata
	*L = (struct pixel *) malloc((*latime_img) * (*inaltime_img) * sizeof(struct pixel));

	//Citirea imaginii in forma liniarizata
	if (okType)
		printf("\nCitire...");
	unsigned char *tmp = (unsigned char*)malloc(3);
	for (i = 0; i < (*inaltime_img); i++)
	{
		for (j = 0; j < (*latime_img); j++)
		{
			fread(&(*L)[((*inaltime_img) - 1 - i)*(*latime_img) + j].B, 1, 1, fin);
			fread(&(*L)[((*inaltime_img) - 1 - i)*(*latime_img) + j].G, 1, 1, fin);
			fread(&(*L)[((*inaltime_img) - 1 - i)*(*latime_img) + j].R, 1, 1, fin);
		}
		if (*padding)
			fread(tmp, 1, (*padding), fin);
		if (i % 128 == 0 && okType) //vizualizare progres
			printf(".");
	}
	printf("\n");
	fclose(fin);
	free(tmp);
	return dim_img;
}

void saveImageLin_C3(char *nume_fisier_destinatie, unsigned char *header, struct pixel *L, unsigned int latime_img, unsigned int inaltime_img, unsigned int padding)
{
	FILE *fout;
	fout = fopen(nume_fisier_destinatie, "wb");
	int i, j;

	//Salvarea header ului imaginii
	for (i = 0; i < 54; i++)
		fwrite(&header[i], 1, 1, fout);

	//Salvarea imaginii din forma liniarizata
	unsigned char *tmp = (unsigned char*)malloc(3); tmp[0] = 0; tmp[1] = 0; tmp[2] = 0;
	printf("\nSalvare...");
	for (i = 0; i < inaltime_img; i++)
	{
		for (j = 0; j < latime_img; j++)
		{
			fwrite(&L[(inaltime_img - 1 - i)*latime_img + j].B, 1, 1, fout);
			fwrite(&L[(inaltime_img - 1 - i)*latime_img + j].G, 1, 1, fout);
			fwrite(&L[(inaltime_img - 1 - i)*latime_img + j].R, 1, 1, fout);
		}
		if (padding)
			fwrite(tmp, 1, padding, fout);
		if (i % 128 == 0) //vizualizare progres
			printf(".");
	}
	printf("\n");
	fclose(fout);
	free(tmp);
}

void chiSquaredTest_C6(char *nume_fisier_testat)
{
	struct frecventa {
		unsigned int R, G, B;
	};

	unsigned int dim_img, latime_img, inaltime_img, padding;
	unsigned char *header;
	float chiR = 0, chiG = 0, chiB = 0, frecvEstT;
	struct pixel *L;
	int i;

	//Citire imagine testata
	dim_img = loadImageLin_C2(nume_fisier_testat, &header, &L, &latime_img, &inaltime_img, &padding, 0);
	if (!dim_img)
		return;

	FILE *fout = fopen("chi-squared-test-values.txt", "ab");

	struct frecventa *frecv = (struct frecventa*)malloc(256 * sizeof(struct frecventa));
	for (i = 0; i < 256; i++)
	{
		frecv[i].R = 0;
		frecv[i].G = 0;
		frecv[i].B = 0;
	}

	//Gasire frecventa fiecarui canal de culoare pt imaginea testata
	for (i = 0; i < latime_img * inaltime_img; i++)
	{
		frecv[L[i].R].R++;
		frecv[L[i].G].G++;
		frecv[L[i].B].B++;
	}

	//Calculare frecventa estimata teoretic (medie)
	frecvEstT = (latime_img * inaltime_img) / 256.0;

	//Calculare valori chi squared pt fiecare canal de culoare
	for (i = 0; i < 256; i++)
	{
		chiR = chiR + ((frecv[i].R - frecvEstT) * (frecv[i].R - frecvEstT)) / frecvEstT;
		chiG = chiG + ((frecv[i].G - frecvEstT) * (frecv[i].G - frecvEstT)) / frecvEstT;
		chiB = chiB + ((frecv[i].B - frecvEstT) * (frecv[i].B - frecvEstT)) / frecvEstT;
	}

	printf("Chi-squared test on RGB channels for %s:\nR: %.2f\nG: %.2f\nB: %.2f\n", nume_fisier_testat, chiR, chiG, chiB);
	fprintf(fout, "Chi-squared test on RGB channels for %s:\r\nR: %.2f\r\nG: %.2f\r\nB: %.2f\r\n\r\n", nume_fisier_testat, chiR, chiG, chiB);
	fclose(fout);
	free(header);
	free(frecv);
}

void encryptImage_C4(char *nume_fisier_sursa, char *nume_fisier_destinatie, char *secret_key)
{
	unsigned int dim_img, latime_img, inaltime_img, padding, *R, *Perm, R0, SV;
	unsigned char *header;
	struct pixel *L, *LP, *C;
	int i, k;

	//Citire chei secrete
	if (!getSecretKeys(secret_key, &R0, &SV))
		return;

	//Citire imagine initiala
	dim_img = loadImageLin_C2(nume_fisier_sursa, &header, &L, &latime_img, &inaltime_img, &padding, 1);
	if (!dim_img)
		return;

	//Generare secventa de nr aleatoare
	xorshift32_C1(&R, 2 * latime_img * inaltime_img, R0);

	//Generare permutare aleatorie
	pDurstenfeld(&Perm, latime_img*inaltime_img, R);

	//Permutarea pixelilor imaginii
	LP = (struct pixel*) malloc(latime_img * inaltime_img * sizeof(struct pixel));
	for (i = 0; i < latime_img * inaltime_img; i++)
		LP[Perm[i]] = L[i];

	//Schimbarea valorilor pixelilor
	printf("\nCriptare...");
	C = (struct pixel*) malloc(latime_img * inaltime_img * sizeof(struct pixel));
	int sv = SV;

	C[0].B = (sv & 255) ^ LP[0].B ^ (R[latime_img * inaltime_img] & 255);

	sv = sv >> 8;
	R[latime_img * inaltime_img] = R[latime_img * inaltime_img] >> 8;

	C[0].G = (sv & 255) ^ LP[0].G ^ (R[latime_img * inaltime_img] & 255);

	sv = sv >> 8;
	R[latime_img * inaltime_img] = R[latime_img * inaltime_img] >> 8;

	C[0].R = (sv & 255) ^ LP[0].R ^ (R[latime_img * inaltime_img] & 255);

	for (k = 1; k < latime_img * inaltime_img; k++)
	{
		C[k].B = (C[k - 1].B) ^ (LP[k].B) ^ (R[latime_img * inaltime_img + k] & 255);

		R[latime_img * inaltime_img + k] = R[latime_img * inaltime_img + k] >> 8;

		C[k].G = (C[k - 1].G) ^ (LP[k].G) ^ (R[latime_img * inaltime_img + k] & 255);

		R[latime_img * inaltime_img + k] = R[latime_img * inaltime_img + k] >> 8;

		C[k].R = (C[k - 1].R) ^ (LP[k].R) ^ (R[latime_img * inaltime_img + k] & 255);

		if (k % (128 * latime_img) == 0) //vizualizare progres
			printf(".");
	}
	printf("\n");

	//Salvare imagine criptata
	saveImageLin_C3(nume_fisier_destinatie, header, C, latime_img, inaltime_img, padding);

	//Eliberare memorie
	free(L);
	free(LP);
	free(C);
	free(R);
	free(Perm);
	free(header);
	printf("\nImaginea a fost criptata.\n");
	chiSquaredTest_C6(nume_fisier_sursa);
	chiSquaredTest_C6(nume_fisier_destinatie);
}

void decryptImage_C5(char *nume_fisier_decriptat, char *nume_fisier_criptat, char *secret_key)
{
	unsigned int dim_img, latime_img, inaltime_img, padding, *R, *Perm, R0, SV;
	unsigned char *header;
	struct pixel *C, *CP, *L;
	int i;

	//Citire chei secrete
	if (!getSecretKeys(secret_key, &R0, &SV))
		return;

	//Citire imagine initiala
	dim_img = loadImageLin_C2(nume_fisier_criptat, &header, &C, &latime_img, &inaltime_img, &padding, 1);
	if (!dim_img)
		return;

	//Generare secventa de nr aleatoare
	xorshift32_C1(&R, 2 * latime_img * inaltime_img, R0);

	//Generare permutare aleatorie
	pDurstenfeld(&Perm, latime_img * inaltime_img, R);

	//Calculare inversa permutarii
	inversaPerm(&Perm, latime_img * inaltime_img);

	//Restaurarea valorilor pixelilor
	printf("\nDecriptare...");
	CP = (struct pixel*) malloc(latime_img * inaltime_img * sizeof(struct pixel));
	int sv = SV;

	CP[0].B = (sv & 255) ^ (C[0].B) ^ (R[latime_img * inaltime_img] & 255);

	sv = sv >> 8;
	R[latime_img * inaltime_img] = R[latime_img * inaltime_img] >> 8;

	CP[0].G = (sv & 255) ^ (C[0].G) ^ (R[latime_img * inaltime_img] & 255);

	sv = sv >> 8;
	R[latime_img * inaltime_img] = R[latime_img * inaltime_img] >> 8;

	CP[0].R = (sv & 255) ^ (C[0].R) ^ (R[latime_img * inaltime_img] & 255);

	for (int k = 1; k < latime_img * inaltime_img; k++)
	{
		CP[k].B = (C[k - 1].B) ^ (C[k].B) ^ (R[latime_img * inaltime_img + k] & 255);

		R[latime_img * inaltime_img + k] = R[latime_img * inaltime_img + k] >> 8;

		CP[k].G = (C[k - 1].G) ^ (C[k].G) ^ (R[latime_img * inaltime_img + k] & 255);

		R[latime_img * inaltime_img + k] = R[latime_img * inaltime_img + k] >> 8;

		CP[k].R = (C[k - 1].R) ^ (C[k].R) ^ (R[latime_img * inaltime_img + k] & 255);

		if (k % (128 * latime_img) == 0) //vizualizare progres
			printf(".");
	}
	printf("\n");

	//Permutarea pixelilor imaginii
	L = (struct pixel*) malloc(latime_img * inaltime_img * sizeof(struct pixel));
	for (i = 0; i < latime_img * inaltime_img; i++)
		L[Perm[i]] = CP[i];

	//Salvare imagine decriptata
	saveImageLin_C3(nume_fisier_decriptat, header, L, latime_img, inaltime_img, padding);

	//Eliberare memorie
	free(L);
	free(CP);
	free(C);
	free(R);
	free(Perm);
	free(header);
	printf("\nImaginea a fost decriptata.\n");
}

void setColors(struct pixel **C)
{
	(*C) = (struct pixel*) malloc(10 * sizeof(struct pixel));
	(*C)[0].R = 255; (*C)[0].G = 0;   (*C)[0].B = 0;
	(*C)[1].R = 255; (*C)[1].G = 255; (*C)[1].B = 0;
	(*C)[2].R = 0;   (*C)[2].G = 255; (*C)[2].B = 0;
	(*C)[3].R = 0;   (*C)[3].G = 255; (*C)[3].B = 255;
	(*C)[4].R = 255; (*C)[4].G = 0;   (*C)[4].B = 255;
	(*C)[5].R = 0;   (*C)[5].G = 0;   (*C)[5].B = 255;
	(*C)[6].R = 192; (*C)[6].G = 192; (*C)[6].B = 192;
	(*C)[7].R = 255; (*C)[7].G = 140; (*C)[7].B = 0;
	(*C)[8].R = 128; (*C)[8].G = 0;   (*C)[8].B = 128;
	(*C)[9].R = 128; (*C)[9].G = 0;   (*C)[9].B = 0;
}

struct pixel **loadMatI(char *nume_fisier_sursa, unsigned char **header)
{
	unsigned int latime_img, inaltime_img, padding;
	struct pixel **MatI;
	int i, j;

	//Deschidere fisier imagine originala
	FILE *f = fopen(nume_fisier_sursa, "rb"); //nu poate fi NULL

	//Latime si inaltime imagine
	fseek(f, 18, SEEK_SET);
	fread(&latime_img, sizeof(unsigned int), 1, f);
	fread(&inaltime_img, sizeof(unsigned int), 1, f);

	//Calculare padding
	if (latime_img % 4 != 0)
		padding = 4 - (3 * latime_img) % 4;
	else
		padding = 0;

	//Copierea header-ului
	*header = (unsigned char*)malloc(54);
	fseek(f, 0, SEEK_SET);
	for (i = 0; i < 54; i++)
		fread(&(*header)[i], 1, 1, f);

	//Citire imagine
	unsigned char *tmp = (unsigned char*)malloc(3);
	printf("\nIncarcare..");
	MatI = (struct pixel **) malloc(inaltime_img * sizeof(struct pixel*));
	for (i = 0; i < inaltime_img; i++)
	{
		MatI[i] = (struct pixel *) malloc(latime_img * sizeof(struct pixel));
		for (j = 0; j < latime_img; j++)
		{
			fread(&MatI[i][j].B, 1, 1, f);
			fread(&MatI[i][j].G, 1, 1, f);
			fread(&MatI[i][j].R, 1, 1, f);
		}
		if (padding)
			fread(tmp, 1, padding, f);
		if (i % 256 == 0) //vizualizare progres
			printf(".");
	}
	printf(" si ");
	fclose(f);
	free(tmp);
	return MatI;
}

void saveMatI(char *nume_fisier_dest, struct pixel **MatI, unsigned char *header)
{
	unsigned int latime_img, inaltime_img, padding;
	int i, j;

	FILE *fout = fopen(nume_fisier_dest, "wb");

	//Latime si inaltime imagine
	latime_img = *(unsigned int*)&header[18];
	inaltime_img = *(unsigned int*)&header[22];

	//Calculare padding
	if (latime_img % 4 != 0)
		padding = 4 - (3 * latime_img) % 4;
	else
		padding = 0;

	//Salvare imagine cu pattern uri recunoscute
	fwrite(header, 1, 54, fout);
	unsigned char *tmp = (unsigned char*)malloc(3); tmp[0] = 0; tmp[1] = 0; tmp[2] = 0;
	printf("salvare..");
	for (i = 0; i < inaltime_img; i++)
	{
		for (j = 0; j < latime_img; j++)
		{
			fwrite(&MatI[i][j].B, 1, 1, fout);
			fwrite(&MatI[i][j].G, 1, 1, fout);
			fwrite(&MatI[i][j].R, 1, 1, fout);
		}
		if (padding)
			fwrite(tmp, 1, padding, fout);
		if (i % 256 == 0) //vizualizare progres
			printf(".");
	}
	printf("\n");
	free(tmp);
	fclose(fout);
}

struct detectie windowMatch(struct pixel **MatI, struct pixel **S, int x, int y, unsigned int inaltime, unsigned int latime, int nrSab)
{
	struct detectie det;
	int i, j;
	float n, medS, dS, medF, dF;

	//det.corr = ?
	det.sablon = nrSab;
	det.inaltime = inaltime;
	det.latime = latime;
	det.x = x;
	det.y = y;

	n = latime * inaltime;
	//Calculare media pixelilor in sablon (medS)
	medS = 0;
	for (i = 0; i < inaltime; i++)
		for (j = 0; j < latime; j++)
			medS += S[i][j].R;
	medS = medS / n;

	//Calculare media pixelilor in fereastra imaginii (medF)
	medF = 0;
	for (i = x; i < x + inaltime; i++)
		for (j = y; j < y + latime; j++)
			medF += MatI[i][j].R;
	medF = medF / n;

	//Calculare deviatie standard in sablon (dS)
	dS = 0;
	for (i = 0; i < inaltime; i++)
		for (j = 0; j < latime; j++)
			dS += (S[i][j].R - medS) * (S[i][j].R - medS);
	dS = sqrt((1 / (n - 1)) * dS);

	//Calculare deviatie standard in fereastra imaginii (dF)
	dF = 0;
	for (i = x; i < x + inaltime; i++)
		for (j = y; j < y + latime; j++)
			dF += (MatI[i][j].R - medF) * (MatI[i][j].R - medF);
	dF = sqrt((1 / (n - 1)) * dF);

	//Calculare corelatie (det.corr)
	det.corr = 0;
	for (i = 0; i < inaltime; i++)
		for (j = 0; j < latime; j++)
			det.corr += (1 / (dF * dS)) * (MatI[i + x][j + y].R - medF) * (S[i][j].R - medS);
	det.corr = (1 / n) * det.corr;

	return det;
}

struct detectie *templateMatching_C7(char *nume_fisier_sursa, char *nume_sablon, float prag, unsigned int *cnt, int nrSab)
{
	unsigned int latime_img, inaltime_img, padding, latime_sab, inaltime_sab, padding_sab;
	char *nume_sablon_grayscale;
	struct detectie *arrDet;
	struct pixel **MatI, **S;
	int i, j;
	FILE *f, *fS;

	//Deschidere fisier imagine
	f = fopen(nume_fisier_sursa, "rb");
	if (f == NULL)
	{
		printf("Nu am gasit imaginea sursa din care citesc.\n");
		return NULL;
	}

	//Latime si inaltime imagine
	fseek(f, 18, SEEK_SET);
	fread(&latime_img, sizeof(unsigned int), 1, f);
	fread(&inaltime_img, sizeof(unsigned int), 1, f);

	//Calculare padding
	if (latime_img % 4 != 0)
		padding = 4 - (3 * latime_img) % 4;
	else
		padding = 0;

	//Citire imagine
	fseek(f, 54, SEEK_SET);
	unsigned char *tmp = (unsigned char*)malloc(3);
	MatI = (struct pixel **) malloc(inaltime_img * sizeof(struct pixel*));
	for (i = 0; i < inaltime_img; i++)
	{
		MatI[i] = (struct pixel *) malloc(latime_img * sizeof(struct pixel));
		for (j = 0; j < latime_img; j++)
		{
			fread(&MatI[i][j].B, 1, 1, f);
			fread(&MatI[i][j].G, 1, 1, f);
			fread(&MatI[i][j].R, 1, 1, f);
		}
		if (padding)
			fread(tmp, 1, padding, f);
	}
	fclose(f);

	//Transformare grayscale a sablonului
	nume_sablon_grayscale = (char*)malloc(60);
	strcpy(nume_sablon_grayscale, "grayscale_");
	strcat(nume_sablon_grayscale, nume_sablon);
	grayscaleImage(nume_sablon, nume_sablon_grayscale, 0);

	//Deschidere fisier sablon
	fS = fopen(nume_sablon_grayscale, "rb");
	if (fS == NULL)
	{
		printf("Nu am gasit imaginea sablonului din care citesc.\n");
		return NULL;
	}

	//Latime si inaltime sablon
	fseek(fS, 18, SEEK_SET);
	fread(&latime_sab, sizeof(unsigned int), 1, fS);
	fread(&inaltime_sab, sizeof(unsigned int), 1, fS);

	//Calculare padding
	if (latime_sab % 4 != 0)
		padding_sab = 4 - (3 * latime_sab) % 4;
	else
		padding_sab = 0;

	//Citire sablon grayscale
	fseek(fS, 54, SEEK_SET);
	S = (struct pixel **) malloc(inaltime_sab * sizeof(struct pixel*));
	for (i = 0; i < inaltime_sab; i++)
	{
		S[i] = (struct pixel *) malloc(latime_sab * sizeof(struct pixel));
		for (j = 0; j < latime_sab; j++)
		{
			fread(&S[i][j].B, 1, 1, fS);
			fread(&S[i][j].G, 1, 1, fS);
			fread(&S[i][j].R, 1, 1, fS);
		}
		if (padding_sab)
			fread(tmp, 1, padding_sab, fS);
	}
	fclose(fS);

	//Glisarea sablonului si retinerea corelatiile peste pragul indicat (0.5)
	arrDet = (struct detectie*) malloc(inaltime_img * latime_img * sizeof(struct detectie));
	*cnt = 0;
	for (i = 0; i <= inaltime_img - inaltime_sab; i++)
	{
		for (j = 0; j <= latime_img - latime_sab; j++)
		{
			arrDet[*cnt] = windowMatch(MatI, S, i, j, inaltime_sab, latime_sab, nrSab);
			if (arrDet[*cnt].corr >= prag)
				(*cnt)++;
		}
	}

	//Eliberare memorie
	for (i = 0; i < inaltime_img; i++)
		free(MatI[i]);
	free(MatI);
	for (i = 0; i < inaltime_sab; i++)
		free(S[i]);
	free(S);
	free(nume_sablon_grayscale);

	return arrDet;
}

void drawBorder_C8(struct pixel **MatI, struct detectie fI, struct pixel C)
{
	int i;
	for (i = fI.x; i < fI.x + fI.inaltime; i++) //stanga si dreapta
	{
		MatI[i][fI.y].B = C.B;
		MatI[i][fI.y].G = C.G;
		MatI[i][fI.y].R = C.R;
		MatI[i][fI.y + fI.latime - 1].B = C.B;
		MatI[i][fI.y + fI.latime - 1].G = C.G;
		MatI[i][fI.y + fI.latime - 1].R = C.R;
	}
	for (i = fI.y; i < fI.y + fI.latime; i++) //sus si jos
	{
		MatI[fI.x][i].B = C.B;
		MatI[fI.x][i].G = C.G;
		MatI[fI.x][i].R = C.R;
		MatI[fI.x + fI.inaltime - 1][i].B = C.B;
		MatI[fI.x + fI.inaltime - 1][i].G = C.G;
		MatI[fI.x + fI.inaltime - 1][i].R = C.R;
	}
}

int cmp_C9(const void *a, const void *b)
{
	if ((*(struct detectie*)b).corr - (*(struct detectie*)a).corr > 0.00)
		return 1;
	else
		return -1;
}

void sortDetect_C9(struct detectie *arr, unsigned int n)
{
	qsort(arr, n, sizeof(struct detectie), cmp_C9);
}

float overlap(struct detectie a, struct detectie b)
{
	float Aa, Ab, Aanb, Aaub;
	if (a.x + a.inaltime <= b.x || a.y + a.latime <= b.y) //nu se intersecteaza
		return 0;
	if (b.x + b.inaltime <= a.x || b.y + b.latime <= a.y) //nu se intersecteaza
		return 0;
	Aa = a.inaltime * a.latime;  //Arie detectie a
	Ab = b.inaltime * b.latime;  //Arie detectie b
	Aanb = (min(a.inaltime, b.inaltime) - abs(a.x - b.x)) * (min(a.latime, b.latime) - abs(a.y - b.y)); //Arie intersectie a n b
	Aaub = Aa + Ab - Aanb;   //Arie reuniune a u b
	return Aanb / Aaub;      //Factor de suprapunere
}

void delete_non_max_C10(struct detectie *arr, unsigned int *n)
{
	int i, j, nAux = 0;
	struct detectie *arrAux = (struct detectie*)malloc((*n) * sizeof(struct detectie));

	for (i = 0; i < (*n); i++)
	{
		if (arr[i].corr == 0) //daca a fost sters
			continue;
		arrAux[nAux] = arr[i];
		nAux++;
		for (j = i + 1; j < (*n); j++)
			if (overlap(arr[i], arr[j]) > 0.2) //verificare suprapunere
				arr[j].corr = 0; //eliminare
	}
	for (i = 0; i < nAux; i++) //actualizare vector cu detectii
		arr[i] = arrAux[i];
	(*n) = nAux; //actualizare nr detectii
	free(arrAux);
}

void patternMatching(char *nume_fisier_sursa, char *nume_fisier_rec, char *nume_fisier_sabloane, float ps)
{
	char *nume_fisier_grayscale, **sabloane;
	unsigned int dimImg, cntDetFi, cntDet;
	unsigned char *header;
	struct detectie *arrFiSab, *arrDetected;
	struct pixel **MatI, *culoare;
	int i, j;

	//Deschidere fisier cu adresele sabloanelor
	FILE *fSab = fopen(nume_fisier_sabloane, "rb");
	if (fSab == NULL)
	{
		printf("Fisierul pt sabloane nu a fost gasit.\n");
		return;
	}

	//Preluarea adreselor sabloanelor
	sabloane = (char**)malloc(10 * sizeof(char*));
	for (i = 0; i < 10; i++)
	{
		sabloane[i] = (char*)malloc(65);
		fscanf(fSab, "%s", sabloane[i]);
	}
	fclose(fSab);

	nume_fisier_grayscale = (char*)malloc(100);
	strcpy(nume_fisier_grayscale, "grayscale_");
	strcat(nume_fisier_grayscale, nume_fisier_sursa);

	//Transformare grayscale a imaginii sursa
	dimImg = grayscaleImage(nume_fisier_sursa, nume_fisier_grayscale, 1);
	if (!dimImg)
	{
		free(nume_fisier_grayscale);
		return;
	}

	//Gasire detectii pt fiecare sablon
	arrDetected = (struct detectie*)malloc(dimImg * 10 * sizeof(struct detectie));
	arrFiSab = (struct detectie*)malloc(dimImg * sizeof(struct detectie));
	cntDet = 0;
	printf("\nCautare detectii...");
	for (i = 0; i < 10; i++)
	{
		//cntDetFi = 0;
		arrFiSab = templateMatching_C7(nume_fisier_grayscale, sabloane[i], ps, &cntDetFi, i);
		//Adaugare a detectiilor gasite pt fiecare sablon in parte la un loc
		for (j = 0; j < cntDetFi; j++)
		{
			arrDetected[cntDet] = arrFiSab[j];
			cntDet++;
		}
		//if (i % 2)
		printf(".");
	}
	printf("\n");
	//Eliberare memorie #1
	for (i = 0; i < 10; i++)
		free(sabloane[i]);
	free(sabloane);
	free(arrFiSab);
	free(nume_fisier_grayscale);

	//Sortare descrescatoare a detectiilor
	sortDetect_C9(arrDetected, cntDet);

	//Eliminarea non-maximelor
	delete_non_max_C10(arrDetected, &cntDet);

	//Citire imagine originala in matrice
	MatI = loadMatI(nume_fisier_sursa, &header);

	//Setare culori pt desenare contur
	setColors(&culoare);

	//Desenare bordura pt fiecare cifra recunoscuta
	for (i = 0; i < cntDet; i++)
		drawBorder_C8(MatI, arrDetected[i], culoare[arrDetected[i].sablon]);

	//Salvare imagine cu recunoasterea finalizara
	saveMatI(nume_fisier_rec, MatI, header);

	printf("\nRecunoasterea de pattern-uri a fost efectuata.\n");

	//Eliberare memorie #2
	for (i = 0; i < *(unsigned int*)&header[22]; i++)
		free(MatI[i]);
	free(MatI);
	free(header);
	free(culoare);
	free(arrDetected);
}

void menu_C11()
{
	int t;
	char *nume_img_sursa = (char*)malloc(90);
	while (1)
	{
		printf("1. Criptare imagine\n");
		printf("2. Decriptare imagine\n");
		printf("3. Recunoastere pattern-uri\n");
		printf("0. Exit\n");
		scanf("%d", &t);
		if (!t)
			break;
		printf("Nume fisier imagine: ");
		scanf("%s", nume_img_sursa);
		switch (t)
		{
		case 1:
		{
			char *nume_img_encrypted = (char*)malloc(100);
			strcpy(nume_img_encrypted, "enc_");
			strcat(nume_img_encrypted, nume_img_sursa);
			printf("Nume fisier cheie secreta: ");
			char *nume_cheie_secreta = (char*)malloc(55);
			scanf("%s", nume_cheie_secreta);
			encryptImage_C4(nume_img_sursa, nume_img_encrypted, nume_cheie_secreta);
			free(nume_img_encrypted);
			free(nume_cheie_secreta);
		}
		break;
		case 2:
		{
			char *nume_img_decrypted = (char*)malloc(100);
			strcpy(nume_img_decrypted, "dec_");
			strcat(nume_img_decrypted, nume_img_sursa);
			printf("Nume fisier cheie secreta: ");
			char *nume_cheie_secreta = (char*)malloc(55);
			scanf("%s", nume_cheie_secreta);
			decryptImage_C5(nume_img_decrypted, nume_img_sursa, nume_cheie_secreta);
			free(nume_img_decrypted);
			free(nume_cheie_secreta);
		}
		break;
		case 3:
		{
			float prag = 0.5;
			char *nume_img_rec = (char*)malloc(100);
			strcpy(nume_img_rec, "rec_");
			strcat(nume_img_rec, nume_img_sursa);
			printf("Nume fisier sabloane: ");
			char *nume_fisier_sabloane = (char*)malloc(50);
			scanf("%s", nume_fisier_sabloane);
			//printf("Prag: ");
			//scanf("%f", &prag);
			patternMatching(nume_img_sursa, nume_img_rec, nume_fisier_sabloane, prag);
			free(nume_img_rec);
			free(nume_fisier_sabloane);
		}
		break;
		case 0:
			break;
		default: printf("Optiune invalida.\n");
			break;
		}
		printf("\n");
	}
	free(nume_img_sursa);
}

int main()
{
	menu_C11();
	return 0;
}