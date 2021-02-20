//dvofazni simpleks
//Jovana Radjenovic 93/2015

#include <stdio.h>
#include <stdlib.h>

double* alociraj_i_ucitaj_c(FILE*, int, double*);
double** alociraj_i_ucitaj_A(FILE*, int, int, double**, double*);
double* alociraj_i_ucitaj_b(FILE*, int, double*, double*, double**);
int* alociraj_baza(int*, int, int, double*, double**, double*);
int* alociraj_kec(int*, int, int, double*, double**, double*, int*);
int* alociraj_w(int*, int, int, double*, double**, double*, int*, int*);
double** alociraj_Awb(double**, int,  double*, double**, double*, int*, int*, int*);
double* alociraj_Awbi(double**, int, int, int, int, double*, double**, double*, int*, int*, int*);
double* alociraj_tab(double*, int, int, int, double*, double**, double*, int*, int*, int*, double**);
int* alociraj_izbaceni(int*, int, int, double*, double**, double*, int*, int*, int*, double**, double*);
double** alociraj_Ab(double**, int, int, double*, double**, double*, int*, int*, int*, double**, double*, int*);
double* alociraj_Abi(double**, int, int, int, double*, double**, double*, int*, int*, int*, double**, double*, int*);
int nadji_negativan(int, double*);
void stampaj_tablicu1(FILE*, int, int, int, double*, double**);
void stampaj_tablicu2(FILE*, int, int, double*, double**);
void oslobodi_alocirano(int, double*, double**, double*, int*, int*, int*, double**, double*, int*, double**);

int main() {
	int i, j, p, t;
	int cd, m, n, mm = 0;
	int k, v;
	double min;
	double mnozilac;
	double* c;
	double* b;
	double** A;

	//Ucitavamo podatke

    //Kako uneti podatke u fajl ulaz.txt:
    /* Ako je dat problem oblika (kanonski oblik se navodi; ukoliko problem nije u kanonskom obliku najpre ga svesti na isti a zatim gledati parametre):


		(min) 2x1 + 3x3 + x4
		p.o	-x2 - x3 + x4 = 3
			2x1 + 2x3 + 4x4 = 12
			x1 + x2 + 2x3 + x4 = 3
			x1,...,x4 >= 0

	Podatke treba navesti u sledecem obliku:
		4		//Broj prromenljivih
		2 0 3 1		//c
		3 4		//Dimenzija matrice A
		0 -1 -1 1	//Matrica A
		2 0 2 4
		1 1 2 1
		3 12 3		//b				*/


	//Otvaramo fajl za citanje
	FILE* in = fopen("ulaz.txt", "r");

	//cd - broj promenjivih
	fscanf(in, "%d", &cd);

	c = alociraj_i_ucitaj_c(in, cd, c);

    //Dimenzija matrice A
    fscanf(in, "%d %d", &m, &n);

	A = alociraj_i_ucitaj_A(in, m, n, A, c);
	b = alociraj_i_ucitaj_b(in, m, b, c, A);
	
	//Zatvaramo ulaz
	fclose(in);

	//Zavrseno ucitavanje


    //baza - ima vrednost 1 ako je promenljiva u bazi; 0 ako nije
    //brw - broj vestackih promenljivih
	int* baza = alociraj_baza(baza, n, m, c, A, b);

    int brw = 0;

    //kec - cuva po jednu vrednost za svaku promenljivu
    //1) 0 ako nije u bazi
    //2) ako promenljiva jeste u bazi...vrstu u kojoj se nalazi kec i plus jedan
    // (plus jedan je dodato u slucaju da je kec u nultoj vrsti da ne bi delovalo kao da promenljiva nije u bazi)
    //Alokacija uz proveru
	int* kec = alociraj_kec(kec, n, m, c, A, b, baza);

    //punimo te nizove
    for (i = 0; i < n; i++) {
        baza[i] = 0;
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            if (A[j][i] != 0 && A[j][i] != 1) {
                baza[i] = 0;
                kec[i] = 0;
                break;
            }
            if (A[j][i] == 1 && baza[i] != 1) {
                baza[i] = 1;
                kec[i] = j + 1;
            }
            else if (A[j][i] == 1 && baza[i] == 1) {
                baza[i] = 0;
                kec[i] = 0;
                break;
            }
        }
    }

    for (i = 0; i < n; i++) {
        if (baza[i] == 1)
            brw++;
    }

    //brw je do ovde jednak broju promenljivih koje su u bazi
    //sada postaje broj vestackih promenljivih koje su nam potrebne
    brw = m - brw;

    //w - isto kao kec samo sto je za vestacke promenljive
	int* w = alociraj_w(w, brw, m, c, A, b, baza, kec);

    //punimo w
    for (j = 0; j < brw; j++) {
        for (i = 0; i < n; i++) {
            if (kec[i] == 0) {
                w[j] = kec[i];
                break;
            }
        }
    }

    k = 0;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            if (kec[j] == i + 1) {
                break;
            }
        }
        if (j == n && kec[n - 1] != i + 1) {
            w[k] = i + 1;
            k++;
        }
    }

    //Awb matrica koja odgovara matrici tablice...ona sadrzi matricu A, matricu sa vestackim promenljivim i matricu b
	double** Awb;

    //Alokacija uz proveru
	Awb = alociraj_Awb(Awb, m, c, A, b, baza,  kec, w);

    //punimo Awb
    for (i = 0; i < m; i++) {
    	//Alokacija za vrstu i
	Awb[i] = alociraj_Awbi(Awb, i, brw, m, n, c, A, b, baza, kec, w);

        for (j = 0; j < n; j++) {
            Awb[i][j] = A[i][j];
        }
        for (j = n; j < n + brw; j++) {
            if (i + 1 == w[j - n])
                Awb[i][j] = 1;
            else
                Awb[i][j] = 0;
        }
        Awb[i][n + brw] = b[i];
    }

    // tab cuva vrednosti koeficijenata f-je ciji se min trazi
    //dakle to je onaj deo tabele ispod isprekidane linije (------)
    //Alokacija uz proveru
	double* tab = alociraj_tab(tab, brw, n, m, c, A, b, baza, kec, w, Awb);

    //punimo tab
    for (i = 0; i < n; i++) {
        tab[i] = 0;
    }

    for (i = n; i < n + brw; i++) {
        tab[i] = 1;
    }

    tab[n + brw] = 0;

	//Otvaramo fajl za ispis
	FILE* out = fopen("izlaz.txt", "w");

    //stampamo tablicu
	stampaj_tablicu1(out, brw, n, m, tab, Awb);


    //na pocetku u tab imamo sve nule osim jedinica koje idu uz vestacke promenjlive
    //vrsimo transformacije nad tablicom tako da dobijemo nule u tab uz vestacke promenjive
    k = 0;
    p = 0;
    int l;
    for (i = 0; i < m; i++) {
        if (p < n && kec[p] != 0) {
            l = kec[p] - 1;
            p++;
            for (j = 0; j < m; j++) {
                for (t = 0; t < n + brw + 1; t++) {
                    double r = Awb[l][t];
                    tab[t] = tab[t] - r;
                }
            }
        }
        if (k < brw && w[k] != 0) {
            l = w[k] - 1;
            k++;
            for (t = 0; t < n + brw + 1; t++) {
                double  r = Awb[l][t];
                tab[t] = tab[t] - r;
            }
        }
    }

    //stampamo tablicu
	stampaj_tablicu1(out, brw, n, m, tab, Awb);

    //izbaceni_iz_baze prati koje promenljive bivaju izbacene iz baze
    //Alokacija uz proveru
	int* izbaceni_iz_baze = alociraj_izbaceni(izbaceni_iz_baze, brw, m, c, A, b, baza, kec, w, Awb, tab) ;

    for (i = 0; i < brw; i++) {
        izbaceni_iz_baze[i] = 0;
    }

    //while petlja se odvija sve dok ima negativnih elemenata u tab
    //najpre se nalazi negativan element sa najmanjim indeksom (ko sto nalaze Blendovo pravilo) i cuva se kolona u kojoj se nalazi (to je k)
    //zatim se trazi min{b[i]/Awb[i][k] | Awb[i][k] > 0} = b[v]/Awb[v][k]...dakle v cuva vrstu koja odgovara min
    //obratiti paznju...u matricu Awb je ukljucena i matrica b kao poslednja kolona medjutim f-ji nadji_negativan je kao broj kolona prosledjeno
    //n+brw a ne n+brw+1, pa to nije problem
    //nasli smo dakle pivot on je u vrsti v i koloni k...tu vrstu delimo sa vrednoscu pivota...pivot postaje 1
    //zatim vrsimo transformacije tablice kako bismo dobili nule ispod i iznad pivota...tj. ubacujemo promenljivu u bazu
    while ((k = nadji_negativan(n + brw, tab)) != -1) {
        min = -1;

        //trazimo min
        for (i = 0; i < m; i++) {
            if (Awb[i][k] > 0) {
                min = b[i] * 1.0 / Awb[i][k];
                v = i;
                break;
            }
        }

        for (j = v + 1; j < m; j++) {
            if (Awb[j][k] > 0) {
                if (b[j] * 1.0 / Awb[j][k] < min) {
                    min = b[j] * 1.0 / Awb[j][k];
                    v = j;
                }
            }
        }

        //pratimo da li neka od vestackih promenljivih izlazi iz baze
        //ako ne onda neka promenljiva izlazi iz baze...to je sledeci for
        for (i = 0; i < brw; i++) {
            if (w[i] == v + 1) {
                izbaceni_iz_baze[i] = 1;
                mm++;
            }
        }


        //promenljiva koja je u bazi a kec joj je u vrsti koja odgovara pivotu biva izbacena iz baze
        for (i = 0; i < n; i++) {
            if (kec[i] == v + 1) {
                kec[i] = 0;
                break;
            }
        }

        //promenljiva kojoj odgovara pivot biva ubacena u bazu
        kec[k] = v + 1;

        //ako nema pozitivnih elemenata u vrsti v problem je neogranicen
        if (min < 0) {
            fprintf(out, "Problem je neogranicen!");
            return 0;
        }

        //stampamo min na dve decimale
        fprintf(out, "\nmin %.2f\n", min);

        //delimo vrstu u kojoj je pivot
        double s = Awb[v][k];
        for (j = 0; j < n + brw + 1; j++) {
            Awb[v][j] = Awb[v][j] * 1.0 / s;
        }

        //vrsimo transformacije...ispod i iznad pivota dobijamo nule
        for (i = 0; i < m; i++) {
            if (i != v) {
                for (j = 0; j < n + brw + 1; j++) {
                    if (j == 0) {
                        mnozilac = -1 * Awb[i][k];
                    }
                    Awb[i][j] = Awb[i][j] + Awb[v][j] * 1.0 * mnozilac;
                }
            }
        }

        for (j = 0; j < n + brw + 1; j++) {
            if (j == 0) {
                mnozilac = -1 * tab[k];
            }
            tab[j] = tab[j] + Awb[v][j] * 1.0 * mnozilac;
        }

	    //stampamo tablicu
		stampaj_tablicu1(out, brw, n, m, tab, Awb);

    }
    //kraj while petlje

    //mm ce biti broj vrsta u novoj matrici koju cemo formirati
    if (brw > 0) {
        mm = m - (brw - mm);
    }
    else {
        mm = m;
    }

    fprintf(out, "\n\n********************DRUGA FAZA*************************\n\n");

    if (tab[n + brw] > 0.000001) {
        fprintf(out, "\nPocetni problem nema dopustivih tacaka!\n");
        return 0;
    }

    //formiramo novu tablicu u nastavku...tab sada postaje c koje je uneo korisnik
    for (i = 0; i < n; i++) {
        tab[i] = c[i];
    }

    tab[n] = 0;

    int temp[n];
    for (t = 0; t < n; t++) {
        temp[t] = kec[t];
    }

    //vestacke promenljive koje su izbacene iz baze se samo izbace iz tablice
    for (j = n; j < n + brw; j++) {
        if (izbaceni_iz_baze[j - n] == 1) {
            for (i = 0; i < m; i++)
                Awb[i][j] = 0;

            tab[j] = 0;
        }
    }

    for (t = 0; t < n; t++) {
        kec[t] = temp[t];
    }

    //stampamo tablicu
	stampaj_tablicu1(out, brw, n, m, tab, Awb);

    //umesto matrice Awb sada cemo uzimati Ab
	double** Ab;

    int ind = 0;

    //formiramo Ab na adekvatan nacin...izbacivanjem vestackih promenljivih po pravilima izbacivanja
    for (i = 0; i < brw; i++) {
        if (w[i] != 0) {
            int o = w[i];
            v = o - 1;

            if (v < m && izbaceni_iz_baze[i] == 0) {
                for (k = 0; k < n + brw + 1; k++) {
                    double h = Awb[v][k];
                    if (Awb[v][k] > 0.00001 && k != i + n) {
                        ind = 1;
                        for (j = 0; j < n + brw + 1; j++) {
                            Awb[v][j] = Awb[v][j] * 1.0 / h;
                        }

                        for (t = 0; t < m; t++) {
                            if (t != v) {
                                for (p = 0; p < n + brw + 1; p++) {
                                    if (p == 0) {
                                        mnozilac = -1 * Awb[p][k];
                                    }
                                    Awb[t][p] = Awb[t][p] + h * 1.0 * mnozilac;
                                }
                            }
                        }

                        izbaceni_iz_baze[i] = 1;
                        for (t = v; t < n; t++) {
                            kec[t]--;
                        }
                        break;
                    }
                }

                if (ind == 0 && k == n + brw + 1) {
                    for (t = 0; t < m; t++) {
                        if (kec[t] > v) {
                            kec[t]--;
                        }
                    }
                }
            }
        }
    }


	//Alokacija uz proveru
	Ab = alociraj_Ab(Ab, mm, m, c, A, b, baza, kec, w, Awb, tab, izbaceni_iz_baze);

	//Alokacija po vrstama
	for(i = 0; i < m;i++) {
		Ab[i] = alociraj_Abi(Ab, i, m, n, c, A, b, baza, kec, w, Awb, tab, izbaceni_iz_baze);
	}

    int v2;
    for (v = 0, v2 = 0; v2 < mm; v++, v2++) {

        int indi = 0;

        for (i = 0; i < brw; i++) {
            if (w[i] == v + 1 && izbaceni_iz_baze[i] == 0) {
                indi = 1;
                v2--;
                break;
            }
        }
        if (indi == 1)
            continue;

        for (j = 0; j < n; j++)
            Ab[v2][j] = Awb[v][j];

        Ab[v2][n] = Awb[v][n + brw];

    }

    m = mm;

    //stampamo tablicu
	stampaj_tablicu2(out, n, m, tab, Ab);

    //odavde krece simpleks...najpre namestamo nule u tab uz promenljive koje su u bazi
    for (i = 0; i < n; i++) {
        if (kec[i] != 0) {
            l = kec[i] - 1;
            double po = tab[i];
            for (t = 0; t < n + 1; t++) {
                tab[t] -= po * Ab[l][t] * 1.0;
            }
        }
    }

    //stampamo tablicu
	stampaj_tablicu2(out, n, m, tab, Ab);

    //novo b
    for (i = 0; i < m; i++) {
        b[i] = Ab[i][n];
    }

    //while petlja koja radi isto sto i prethodna samo sto je prilagodjena novoj tablici
    while ((k = nadji_negativan(n, tab)) != -1) {
        min = -1;


        for (i = 0; i < m; i++) {
            if (Ab[i][k] > 0) {
                min = b[i] * 1.0 / Ab[i][k];
                v = i;
                break;
            }
        }

        for (j = v + 1; j < m; j++) {
            if (Ab[j][k] > 0) {
                if (b[j] * 1.0 / Ab[j][k] < min) {
                    min = b[j] * 1.0 / Ab[j][k];
                    v = j;
                }
            }
        }


        for (i = 0; i < n; i++) {
            if (kec[i] == v + 1)
                kec[i] = 0;
        }

        kec[k] = v + 1;

        if (min < 0) {
            fprintf(out, "Problem je neogranicen!");
            return 0;
        }

        fprintf(out, "\nmin %.2f\n", min);

        double s = Ab[v][k];
        for (j = 0; j < n + 1; j++) {
            Ab[v][j] = Ab[v][j] * 1.0 / s;
        }

    	//stampamo tablicu
		stampaj_tablicu2(out, n, m, tab, Ab);

        for (i = 0; i < m; i++) {
            if (i != v) {
                for (j = 0; j < n + 1; j++) {
                    if (j == 0) {
                        mnozilac = -1 * Ab[i][k];
                    }
                    Ab[i][j] = Ab[i][j] + Ab[v][j] * 1.0 * mnozilac;
                }
            }
        }

        for (j = 0; j < n + 1; j++) {
            if (j == 0) {
                mnozilac = -1 * tab[k];
            }
            tab[j] = tab[j] + Ab[v][j] * 1.0 * mnozilac;
        }

    	//stampamo tablicu
		stampaj_tablicu2(out, n, m, tab, Ab);
    }


    //izasli smo iz while petlje jer nema vise negativnih elemenata u tab...dobijamo minf i optimalno resenje
    //i to ispisujemo
    double  minf = -1 * tab[n];
    fprintf(out, "\nminf = %.2f\n", minf);

    fprintf(out, "\nOptimalno resenje: (");
    double opt[n];
    for (i = 0; i < n; i++) {
        if (kec[i] != 0)
            opt[i] = Ab[kec[i] - 1][n];
        else
            opt[i] = 0;
        if (i != n - 1)
            fprintf(out, " %.2f,", opt[i]);
        else
            fprintf(out, " %.2f)\n", opt[i]);
    }

	//Zatvaramo izlaz
	fclose(out);

	//Oslobadjamo alociranu memoriju
	oslobodi_alocirano(m, c, A, b, baza, kec, w, Awb, tab, izbaceni_iz_baze, Ab);

    return 0;
}

//F-ja koja alocira memoriju za c i ucitava njegove vrednosti
double* alociraj_i_ucitaj_c(FILE* in, int cd, double* c) {
	int i;
	
	//Alokacija uz proveru
	c = (double*)malloc(cd * sizeof(double));
	if (c == NULL) {
		fprintf(stderr, "Greska pri alokaciji memorije za %s\n", "c");
		exit(EXIT_FAILURE);
	}
	
	//Ucitavamo vrednosti za c 
	for (i = 0; i < cd; i++) {
		fscanf(in, "%lf", &c[i]);
	}
	
	return c;
}


//F-ja koja alocira memoriju za A i ucitava njene vrednosti
double** alociraj_i_ucitaj_A(FILE* in, int m, int n, double** A, double* c) {
	int i, j;
	
	//Alokacija uz proveru (ukoliko alokacija nije uspesna ispisujemo poruku o gresci,
	//oslobadjamo prethodno alociranu memoriju, u ovom slucaju to je memorija alocirana za c
	//i zaustavljemo program)
	A = (double**)malloc(m * sizeof(double*));
	if (A == NULL) {
		fprintf(stderr, "Greska pri alokaciji memorije za %s\n", "A");
		
		free(c);
		
		exit(EXIT_FAILURE);
	}

    //Ucitavamo matricu A
    for (i = 0;i < m; i++) {
    	//Alokacija za vrstu i
		A[i] = (double*)malloc(n * sizeof(double));
		if (A[i] == NULL) {
			fprintf(stderr, "Greska pri alokaciji memorije za %s\n", "A[i]");
			
			free(c);
			for (j = 0; j < i; j++)
				free(A[j]);
			free(A);
			
			exit(EXIT_FAILURE);
		}
		
        for (j = 0; j < n; j++) {
            fscanf(in, "%lf", &A[i][j]);
        }
    }
	
	return A;
}

//F-ja koja alocira memoriju za b i ucitava njegove vrednosti
double* alociraj_i_ucitaj_b(FILE* in, int m, double* b, double* c, double** A) {
	int i;
	
	//Alokacija uz proveru (ukoliko alokacija nije uspesna ispisujemo poruku o gresci,
	//oslobadjamo prethodno alociranu memoriju  i zaustavljemo program)
	b = (double*)malloc(m * sizeof(double));
	if (b == NULL) {
		fprintf(stderr, "Greska pri alokaciji memorije za %s\n", "b");
		
		free(c);
		for (i = 0; i < m; i++)
			free(A[i]);
		free(A);
		
		exit(EXIT_FAILURE);
	}
	
	//Ucitavamo b
	for (i = 0; i < m; i++) {
	    fscanf(in, "%lf", &b[i]);
	}
	
	return b;
}

//F-ja koja alocira memoriju za "baza"
int* alociraj_baza(int* baza, int n, int m, double* c, double** A, double* b) {
	int i;
	
	//Alokacija uz proveru (ukoliko alokacija nije uspesna ispisujemo poruku o gresci,
	//oslobadjamo prethodno alociranu memoriju  i zaustavljemo program)
	baza = (int*)malloc(n * sizeof(int));
	if (baza == NULL) {
		fprintf(stderr, "Greska pri alokaciji memorije za %s\n", "baza");
		
		free(c);
		for (i = 0; i < m; i++)
			free(A[i]);
		free(A);
		free(b);
		
		exit(EXIT_FAILURE);
	}
	
	return baza;
}

//F-ja koja alocira memoriju za "kec"
int* alociraj_kec(int* kec, int n, int m, double* c, double** A, double* b, int* baza) {
	int i;
	
	//Alokacija uz proveru (ukoliko alokacija nije uspesna ispisujemo poruku o gresci,
	//oslobadjamo prethodno alociranu memoriju  i zaustavljemo program)
	kec = (int*)malloc(n * sizeof(int));
	if (kec == NULL) {
		fprintf(stderr, "Greska pri alokaciji memorije za %s\n", "kec");
		
		free(c);
		for (i = 0; i < m; i++)
			free(A[i]);
		free(A);
		free(b);
		free(baza);
		
		exit(EXIT_FAILURE);
	}
	
	return kec;	
}

//F-ja koja alocira memoriju za w
int* alociraj_w(int* w, int brw, int m, double* c, double** A, double* b, int* baza, int* kec) {
	int i;
	
	//Alokacija uz proveru (ukoliko alokacija nije uspesna ispisujemo poruku o gresci,
	//oslobadjamo prethodno alociranu memoriju  i zaustavljemo program)
	w = (int*)malloc(brw * sizeof(int));
	if (w == NULL) {
		fprintf(stderr, "Greska pri alokaciji memorije za %s\n", "c");
		
		free(c);
		for (i = 0; i < m; i++)
			free(A[i]);
		free(A);
		free(b);
		free(baza);
		free(kec);
		
		exit(EXIT_FAILURE);
	}
	
	return w;
}

//F-ja koja alocira memoriju za Awb
double** alociraj_Awb(double** Awb, int m,  double* c, double** A, double* b, int* baza, int* kec, int* w) {
	int i;
	
	//Alokacija uz proveru (ukoliko alokacija nije uspesna ispisujemo poruku o gresci,
	//oslobadjamo prethodno alociranu memoriju  i zaustavljemo program)
	Awb = (double**)malloc(m * sizeof(double*));
	if (Awb == NULL) {
		fprintf(stderr, "Greska pri alokaciji memorije za %s\n", "Awb");
		
		free(c);
		for (i = 0; i < m; i++)
			free(A[i]);
		free(A);
		free(b);
		free(baza);
		free(kec);
		free(w);
		
		exit(EXIT_FAILURE);
	}
	
	return Awb;
}

//F-ja koja alocira memoriju za Awb[i]
double* alociraj_Awbi(double** Awb, int i, int brw, int m, int n, double* c, double** A, double* b, int* baza, int* kec, int* w) {
	int j;
	
   	//Alokacija za vrstu i
   	//uz proveru (ukoliko alokacija nije uspesna ispisujemo poruku o gresci,
	//oslobadjamo prethodno alociranu memoriju  i zaustavljemo program)
	Awb[i] = (double*)malloc((n + brw + 1) * sizeof(double));
	if (Awb[i] == NULL) {
		fprintf(stderr, "Greska pri alokaciji memorije za %s\n", "Awb[i]");	
				
		free(c);
		for (j = 0; j < m; j++)
			free(A[j]);
		free(A);
		free(b);
		free(baza);
		free(kec);
		free(w);
		for (j = 0; j < i; j++)
			free(Awb[j]);
		free(Awb);
			
		exit(EXIT_FAILURE);
	}
	
	return Awb[i];
}

//F-ja koja alocira memoriju za tab
double* alociraj_tab(double* tab, int brw, int n, int m, double* c, double** A, double* b, int* baza, int* kec, int* w, double** Awb) {
	int i;

	//Alokacija uz proveru (ukoliko alokacija nije uspesna ispisujemo poruku o gresci,
	//oslobadjamo prethodno alociranu memoriju  i zaustavljemo program)
	tab = (double*)malloc((n + brw + 1) * sizeof(double));
	if (tab == NULL) {
		fprintf(stderr, "Greska pri alokaciji memorije za %s\n", "tab");
		
		free(c);
		for (i = 0; i < m; i++)
			free(A[i]);
		free(A);
		free(b);
		free(baza);
		free(kec);
		free(w);
		for (i = 0; i < m; i++)
			free(Awb[i]);
		free(Awb);
		
		exit(EXIT_FAILURE);
	}
		
	return tab;
}

//F-ja koja alocira memoriju za "izbaceni_iz_baze"
int* alociraj_izbaceni(int* izbaceni_iz_baze, int brw, int m, double* c, double** A, double* b, int* baza, int* kec, int* w, double** Awb, double* tab) {
	int i, j;
	
	//Alokacija uz proveru (ukoliko alokacija nije uspesna ispisujemo poruku o gresci,
	//oslobadjamo prethodno alociranu memoriju  i zaustavljemo program)
	izbaceni_iz_baze = (int*)malloc(brw * sizeof(int));
	if (izbaceni_iz_baze == NULL) {
		fprintf(stderr, "Greska pri alokaciji memorije za %s", "izbaceni_iz_baze.");
		
		free(c);
		for (i = 0; i < m; i++)
			free(A[i]);
		free(A);
		free(b);
		free(baza);
		free(kec);
		free(w);
		for (i = 0; i < m; i++)
			free(Awb[i]);
		free(Awb);
		free(tab);
	
		exit(EXIT_FAILURE);
	}
	
	return izbaceni_iz_baze;
}

//F-ja koja alocira memoriju za Ab
double** alociraj_Ab(double** Ab, int mm, int m, double* c, double** A, double* b, int* baza, int* kec, int* w, double** Awb, double* tab, int* izbaceni_iz_baze) {
	int i;
	
	//Alokacija uz proveru (ukoliko alokacija nije uspesna ispisujemo poruku o gresci,
	//oslobadjamo prethodno alociranu memoriju  i zaustavljemo program)
	Ab = (double**)malloc(mm * sizeof(double*));
	if (Ab == NULL) {
		fprintf(stderr, "Greska pri alokaciji memorije za %s", "Ab.");
		
		free(c);
		for (i = 0; i < m; i++)
			free(A[i]);
		free(A);
		free(b);
		free(baza);
		free(kec);
		free(w);
		for (i = 0; i < m; i++)
			free(Awb[i]);
		free(Awb);
		free(tab);
		free(izbaceni_iz_baze);
		
		exit(EXIT_FAILURE);
	}
	
	return Ab;
}

//F-ja koja alocira memoriju za Ab[i]
double* alociraj_Abi(double** Ab, int i, int m, int n, double* c, double** A, double* b, int* baza, int* kec, int* w, double** Awb, double* tab, int* izbaceni_iz_baze) {
	int j;
	
	//Alokacija za vrstu i
	//uz proveru (ukoliko alokacija nije uspesna ispisujemo poruku o gresci,
	//oslobadjamo prethodno alociranu memoriju  i zaustavljemo program)
	Ab[i] = (double*)malloc((n+1) * sizeof(double));
	if (Ab[i] == NULL) {
		fprintf(stderr, "Greska pri alokaciji memorije za %s", "Ab[v2].");
			
		free(c);
		for (j = 0; j < m; j++)
			free(A[j]);
		free(A);
		free(b);
		free(baza);
		free(kec);
		free(w);
		for (j = 0; j < m; j++)
			free(Awb[j]);
		free(Awb);
		free(tab);
		free(izbaceni_iz_baze);
		for(j = 0; j < i; j++)
			free(Ab[j]);
		free(Ab);
						
		exit(EXIT_FAILURE);
	}
		
	return Ab[i];
}

//F-ja koja ce da trazi negativan element u nizu
int nadji_negativan(int h, double* tab) {
    int k, i;

    for (i = 0; i < h; i++) {
        if (tab[i] < -0.00001) {
            k = i;
            return k;
        }
    }

    return -1;
}

//F-ja za ispis tablice za Awb
void stampaj_tablicu1(FILE* out, int brw, int n, int m, double* tab, double** Awb) {
	int i, j;
	
    for (i = 0; i < m; i++) {
        for (j = 0; j < n + brw + 1; j++) {
            fprintf(out, " %.2f ", Awb[i][j]);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "\n------------------------------------------");
    fprintf(out, "\n");

    for (i = 0; i < n + brw + 1; i++) {
        fprintf(out, " %.2f ", tab[i]);
    }
    fprintf(out, "\n========================================");
    fprintf(out, "\n");
}

//F-ja za ispis tablice za Ab
void stampaj_tablicu2(FILE* out, int n, int m, double* tab, double** Ab) {
	int i, j;
	
    for (i = 0; i < m; i++) {
        for (j = 0; j < n + 1; j++) {
            fprintf(out, " %.2f ", Ab[i][j]);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "\n------------------------------------------");
    fprintf(out, "\n");

    for (i = 0; i < n + 1; i++) {
        fprintf(out, " %.2f ", tab[i]);
    }
    fprintf(out, "\n========================================");
    fprintf(out, "\n");
}

//F-ja koja oslobadja alociranu memoriju
void oslobodi_alocirano(int m, double* c, double** A, double* b, int* baza, int* kec, int* w, double** Awb, double* tab, int* izbaceni_iz_baze, double** Ab) {
	int i;

	free(c);
	for (i = 0; i < m; i++)
		free(A[i]);
	free(A);
	free(b);
	free(baza);
	free(kec);
	free(w);
	for (i = 0; i < m; i++)
		free(Awb[i]);
	free(Awb);
	free(tab);
	free(izbaceni_iz_baze);
	for (i = 0; i < m; i++)
		free(Ab[i]);
	free(Ab);
}
