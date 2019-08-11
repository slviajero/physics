/* f2ctmp_lang2hofc.f -- translated by f2c (version of 21 October 1993  13:46:10).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

union {
    struct {
	doublereal xr, xk;
    } _1;
    struct {
	doublereal x, xk;
    } _2;
} rg_;

#define rg_1 (rg_._1)
#define rg_2 (rg_._2)

struct {
    real pi;
} const_;

#define const_1 const_

struct {
    real m, b, v0, db, w, omega, mo, xo, mb, mob, f1, f2;
} par_;

#define par_1 par_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__4 = 4;

/* Main program */ MAIN__(void)
{
    /* Format strings */
    static char fmt_1010[] = "(3e12.4)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3;
    olist o__1;

    /* Builtin functions */
    double atan(doublereal);
    integer f_open(olist *), s_wsle(cilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_wsle(void), s_rsle(cilist *), e_rsle(void);
    double sqrt(doublereal), sinh(doublereal), cosh(doublereal), exp(
	    doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i, j;
    static real q[2406]	/* was [3][802] */;
    extern /* Subroutine */ int s_(integer *, real *, real *, real *, real *, 
	    integer *, real *, real *);
    extern real v_(real *);
    static real x[2406]	/* was [3][802] */, z[3];
    static integer k1, l1, n0, l2;
    static real r2[802];
    static integer ia, nb, nc, na;
    static real sa;
    static integer ng;
    static real dq[2406]	/* was [3][802] */, ds, dt, qi[3], qf[3];
    extern real vf_(real *), vg_(real *);
    static real xf[3], dx[2406]	/* was [3][802] */, xi[3], dw[2406]	/* 
	    was [3][802] */, st, tt, sv, zq[3], r2t[802];
    extern real ran_(void);
    static real sav, dwq[2406]	/* was [3][802] */, svv;
    extern real rgau_(void);
    extern /* Subroutine */ int makran_(real *, real *, integer *, real *);
    static logical cwrite, twrite;

    /* Fortran I/O blocks */
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };
    static cilist io___7 = { 0, 5, 1, 0, 0 };
    static cilist io___9 = { 0, 6, 0, 0, 0 };
    static cilist io___10 = { 0, 5, 0, 0, 0 };
    static cilist io___11 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 5, 0, 0, 0 };
    static cilist io___13 = { 0, 6, 0, 0, 0 };
    static cilist io___14 = { 0, 5, 0, 0, 0 };
    static cilist io___15 = { 0, 6, 0, 0, 0 };
    static cilist io___16 = { 0, 5, 0, 0, 0 };
    static cilist io___17 = { 0, 6, 0, 0, 0 };
    static cilist io___18 = { 0, 5, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };
    static cilist io___20 = { 0, 5, 0, 0, 0 };
    static cilist io___21 = { 0, 6, 0, 0, 0 };
    static cilist io___22 = { 0, 6, 0, 0, 0 };
    static cilist io___23 = { 0, 5, 0, 0, 0 };
    static cilist io___25 = { 0, 6, 0, 0, 0 };
    static cilist io___26 = { 0, 5, 0, 0, 0 };
    static cilist io___28 = { 0, 6, 0, 0, 0 };
    static cilist io___29 = { 0, 5, 0, 0, 0 };
    static cilist io___31 = { 0, 6, 0, 0, 0 };
    static cilist io___32 = { 0, 5, 0, 0, 0 };
    static cilist io___34 = { 0, 6, 0, 0, 0 };
    static cilist io___35 = { 0, 6, 0, 0, 0 };
    static cilist io___36 = { 0, 6, 0, 0, 0 };
    static cilist io___37 = { 0, 6, 0, 0, 0 };
    static cilist io___39 = { 0, 10, 0, 0, 0 };
    static cilist io___40 = { 0, 10, 0, 0, 0 };
    static cilist io___41 = { 0, 10, 0, 0, 0 };
    static cilist io___42 = { 0, 10, 0, 0, 0 };
    static cilist io___68 = { 0, 10, 0, 0, 0 };
    static cilist io___69 = { 0, 10, 0, fmt_1010, 0 };
    static cilist io___70 = { 0, 14, 0, 0, 0 };
    static cilist io___71 = { 0, 12, 0, 0, 0 };


/* -------------------------------------------------------- */
/* Fuehrt Langevin Simulation durch um Konfigurationen zu */
/* gewinnen, die nach der Wirkung von einem Teilchen im */
/* Gausspotential verteilt sind. Berechnet wird direkt die */
/* Streulaenge. */
/* -------------------------------------------------------- */
/* Diese Version sampelt nicht die Anfangs- und Endpunkte */
/* durch den Gaussschen Zufallszahlengenerator, sondern */
/* bezieht Anfangs- und Endpunkt mit in die Langevin Simulation */
/* ein. Die Schleife ueber die Anfangskonfigurationen wird */
/* behalten, jedoch hat sie hier eine andere Bedeutung als */
/* in lang1 */
/* ---------------------------------------------------------- */
/* oszillator version des langevin programms */
/* ---------------------------------------------------------- */
/* diese version implementiert die greensfunktion des */
/* oszillators exakt, dh. eine wirkung mit cosh und sinh */
/* termen wird verwendet */
/* ---------------------------------------------------------- */
/* zusaetzlich wird mit gefalteten potential verglichen */
/* ausgegeben wird nurmehr das verhaeltnis der querschnitte */
/* ----------------------------------------------------------- */
/* weitere modifikation: das effektive potential ist nicht */
/* mehr notwendigerweise das gefaltete potetial, es kann */
/* beliebig gewaehlt werden */
/* ----------------------------------------------------------- */

/* externe Funktionen */




/* Parameter des Projectils: */
/* Parameter des Oszi: */
/* Anfangskonfiguration: */
/* zahl der Anfangskonfigs: */

/* Felder fuer Langevin */


/* Zufallsgenerator */


/* globale variablen */


/* parameter */

    const_1.pi = atan(1.f) * 4.f;
    cwrite = FALSE_;
    twrite = TRUE_;

/* files */

    o__1.oerr = 0;
    o__1.ounit = 10;
    o__1.ofnmlen = 4;
    o__1.ofnm = "test";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    o__1.oerr = 0;
    o__1.ounit = 12;
    o__1.ofnmlen = 5;
    o__1.ofnm = "error";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    o__1.oerr = 0;
    o__1.ounit = 14;
    o__1.ofnmlen = 5;
    o__1.ofnm = "ratio";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);

/* Hauptschleife */

L1:

/*  INITIALIZATION OF THE RANDOM NUMBER GENERATOR */
    rg_1.xr = const_1.pi * 1e11;
    rg_1.xk = 0.;
    for (k1 = 1; k1 <= 20; ++k1) {
	tt = ran_();
/* L10: */
    }

/* Eingabe der Parameter */

    s_wsle(&io___5);
    do_lio(&c__9, &c__1, "Modellparameter :", 17L);
    e_wsle();
    s_wsle(&io___6);
    do_lio(&c__9, &c__1, "Zahl der beta-schritte :", 24L);
    e_wsle();
    i__1 = s_rsle(&io___7);
    if (i__1 != 0) {
	goto L9999;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&nb, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L9999;
    }
    i__1 = e_rsle();
    s_wsle(&io___9);
    do_lio(&c__9, &c__1, "db :", 4L);
    e_wsle();
    s_rsle(&io___10);
    do_lio(&c__4, &c__1, (char *)&par_1.db, (ftnlen)sizeof(real));
    e_rsle();
    s_wsle(&io___11);
    do_lio(&c__9, &c__1, "Weite des Potentials :", 22L);
    e_wsle();
    s_rsle(&io___12);
    do_lio(&c__4, &c__1, (char *)&par_1.b, (ftnlen)sizeof(real));
    e_rsle();
    s_wsle(&io___13);
    do_lio(&c__9, &c__1, "Tiefe des Potentials :", 22L);
    e_wsle();
    s_rsle(&io___14);
    do_lio(&c__4, &c__1, (char *)&par_1.v0, (ftnlen)sizeof(real));
    e_rsle();
    s_wsle(&io___15);
    do_lio(&c__9, &c__1, "Masse des Projectils :", 22L);
    e_wsle();
    s_rsle(&io___16);
    do_lio(&c__4, &c__1, (char *)&par_1.m, (ftnlen)sizeof(real));
    e_rsle();
    s_wsle(&io___17);
    do_lio(&c__9, &c__1, "Oszillatorstaerke :", 19L);
    e_wsle();
    s_rsle(&io___18);
    do_lio(&c__4, &c__1, (char *)&par_1.omega, (ftnlen)sizeof(real));
    e_rsle();
    s_wsle(&io___19);
    do_lio(&c__9, &c__1, "Oszi Masse:", 11L);
    e_wsle();
    s_rsle(&io___20);
    do_lio(&c__4, &c__1, (char *)&par_1.mo, (ftnlen)sizeof(real));
    e_rsle();
    s_wsle(&io___21);
    do_lio(&c__9, &c__1, "Langevin Parameter :", 20L);
    e_wsle();
    s_wsle(&io___22);
    do_lio(&c__9, &c__1, "Zahl der gezogenen Konfigurationen :", 36L);
    e_wsle();
    s_rsle(&io___23);
    do_lio(&c__3, &c__1, (char *)&nc, (ftnlen)sizeof(integer));
    e_rsle();
    s_wsle(&io___25);
    do_lio(&c__9, &c__1, "Gleichgewichtskonfigurationen      :", 36L);
    e_wsle();
    s_rsle(&io___26);
    do_lio(&c__3, &c__1, (char *)&ng, (ftnlen)sizeof(integer));
    e_rsle();
    s_wsle(&io___28);
    do_lio(&c__9, &c__1, "Zahl der ausgelassenen Schritte    :", 36L);
    e_wsle();
    s_rsle(&io___29);
    do_lio(&c__3, &c__1, (char *)&n0, (ftnlen)sizeof(integer));
    e_rsle();
    s_wsle(&io___31);
    do_lio(&c__9, &c__1, "Laenge des Zeitintervalls :", 27L);
    e_wsle();
    s_rsle(&io___32);
    do_lio(&c__4, &c__1, (char *)&dt, (ftnlen)sizeof(real));
    e_rsle();
    s_wsle(&io___34);
    do_lio(&c__9, &c__1, "Simulationszeit :", 17L);
    e_wsle();
    s_wsle(&io___35);
    r__1 = nc * n0 * dt;
    do_lio(&c__4, &c__1, (char *)&r__1, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___36);
    do_lio(&c__9, &c__1, "Korrelationszeit :", 18L);
    e_wsle();
    s_wsle(&io___37);
    r__1 = n0 * dt;
    do_lio(&c__4, &c__1, (char *)&r__1, (ftnlen)sizeof(real));
    e_wsle();
/*      print *,'Zahl der Anfangskonfigurationen :' */
/*      read *,na */
    na = 1;
    s_wsle(&io___39);
    do_lio(&c__9, &c__1, "Lang2hof: Version 1.12.1993", 27L);
    e_wsle();
    s_wsle(&io___40);
    do_lio(&c__3, &c__1, (char *)&na, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&nb, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&nc, (ftnlen)sizeof(integer));
    do_lio(&c__4, &c__1, (char *)&dt, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___41);
    do_lio(&c__4, &c__1, (char *)&par_1.m, (ftnlen)sizeof(real));
    do_lio(&c__4, &c__1, (char *)&par_1.b, (ftnlen)sizeof(real));
    do_lio(&c__4, &c__1, (char *)&par_1.v0, (ftnlen)sizeof(real));
    do_lio(&c__4, &c__1, (char *)&par_1.db, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___42);
    do_lio(&c__4, &c__1, (char *)&par_1.mo, (ftnlen)sizeof(real));
    do_lio(&c__4, &c__1, (char *)&par_1.omega, (ftnlen)sizeof(real));
    e_wsle();

/* berechne bestimmte kombinationen von werten im voraus */

    ds = sqrt(dt * 2.f);
    par_1.xo = sqrt(par_1.mo * par_1.omega);
    par_1.w = par_1.b * 2.f * par_1.b;
    par_1.mob = par_1.mo / par_1.db;
    par_1.mb = par_1.m / par_1.db;
    par_1.f1 = par_1.mo * par_1.omega / sinh(par_1.db * par_1.omega);
    par_1.f2 = cosh(par_1.db * par_1.omega);
    sa = 0.f;
    sav = 0.f;
    i__1 = na;
    for (ia = 1; ia <= i__1; ++ia) {
	sv = 0.f;
	svv = 0.f;
/* ----------------------------------------------------------- */
/* Anfangspunkte */

/*     print *,'Anfangswerte samplen :' */
	for (j = 1; j <= 3; ++j) {
	    qi[j - 1] = rgau_() * par_1.xo;
/* L23: */
	}
	for (j = 1; j <= 3; ++j) {
	    qf[j - 1] = rgau_() * par_1.xo;
/* L24: */
	}
	for (j = 1; j <= 3; ++j) {
	    xi[j - 1] = rgau_() * par_1.b + qi[j - 1];
/* L20: */
	}
	for (j = 1; j <= 3; ++j) {
	    xf[j - 1] = rgau_() * par_1.b + qf[j - 1];
/* L21: */
	}
/* -------------------------------------------------------------- */
/* Anfangswert, straight line zwischen Anfangs und Endpunkt */
/* das kann man sicher besser machen ! */

/*       print *,'Anfangswert fuer Langevin Gleichung machen :' */
	i__2 = nb + 1;
	for (i = 0; i <= i__2; ++i) {
	    for (j = 1; j <= 3; ++j) {
		x[j + i * 3 - 1] = (xf[j - 1] - xi[j - 1]) * i / (real) (nb + 
			1) + xi[j - 1];
		q[j + i * 3 - 1] = (qf[j - 1] - qi[j - 1]) * i / (real) (nb + 
			1) + qi[j - 1];
/* L110: */
	    }
/* L100: */
	}
/* ------------------------------------------------------------- */
/* Langevin Simulation */


/* Gleichgewichtslauf */

/*         print *,'Gleichgewichtslauf :',ia */
	i__2 = ng;
	for (l2 = 1; l2 <= i__2; ++l2) {
/*          print *,l2 */

/* ein schritt */

	    makran_(dw, dwq, &nb, &ds);
	    i__3 = nb + 1;
	    for (i = 0; i <= i__3; ++i) {
		s_(&i, x, q, z, zq, &nb, r2, r2t);
		for (j = 1; j <= 3; ++j) {
		    dx[j + i * 3 - 1] = z[j - 1] * dt + dw[j + i * 3 - 1];
		    dq[j + i * 3 - 1] = zq[j - 1] * dt + dwq[j + i * 3 - 1];
/*                  print *,i,j,z(j)*dt */
/* L192: */
		}
/* L191: */
	    }
	    i__3 = nb + 1;
	    for (i = 0; i <= i__3; ++i) {
		for (j = 1; j <= 3; ++j) {
		    x[j + i * 3 - 1] += dx[j + i * 3 - 1];
		    q[j + i * 3 - 1] += dq[j + i * 3 - 1];
/* L194: */
		}
/* L193: */
	    }

/* ende des schrittes */

/*             do 195 i=1,nb */
/*                do 195 j=1,3 */
/*                   print *,i,j,x(j,i),dx(j,i),dw(j,i) */
/* 195          continue */
/* L190: */
	}

/* eigentliche simulation */

/*       print *,'Simulationsbeginn ',ia */
	i__2 = nc;
	for (l1 = 1; l1 <= i__2; ++l1) {
/*          if (mod(l1,1000).eq.0) print *,l1 */

/* auslassung (n0 schritte) */

	    i__3 = n0;
	    for (l2 = 1; l2 <= i__3; ++l2) {

/* ein schritt */

		makran_(dw, dwq, &nb, &ds);
		i__4 = nb + 1;
		for (i = 0; i <= i__4; ++i) {
		    s_(&i, x, q, z, zq, &nb, r2, r2t);
		    for (j = 1; j <= 3; ++j) {
			dx[j + i * 3 - 1] = z[j - 1] * dt + dw[j + i * 3 - 1];
			dq[j + i * 3 - 1] = zq[j - 1] * dt + dwq[j + i * 3 - 
				1];
/* L230: */
		    }
/* L220: */
		}
		i__4 = nb + 1;
		for (i = 0; i <= i__4; ++i) {
		    for (j = 1; j <= 3; ++j) {
			x[j + i * 3 - 1] += dx[j + i * 3 - 1];
			q[j + i * 3 - 1] += dq[j + i * 3 - 1];
/* L250: */
		    }
/* L240: */
		}

/* ende des schrittes */

/* L210: */
	    }

/* Potentialanteil der Wirkung berechnen */

	    st = 0.f;
	    i__3 = nb;
	    for (i = 1; i <= i__3; ++i) {
		st += (v_(&r2[i]) - vf_(&r2t[i])) * par_1.db;
/* L400: */
	    }
	    st += (v_(&r2[nb + 1]) - vf_(&r2t[nb + 1])) * par_1.db / 2;
	    st += (v_(r2) - vf_(r2t)) * par_1.db / 2;
/*        print *,'Wirkung der ',l1,'.  Konf.:',st */
/*        write(10,*) l1,exp(st) */
/* ------------------------------------------------- */
/* Korrekturfaktor vg --> vf */
/* wird hier angebracht */
/* -------------------------------------------------- */
	    sv += exp(st) * (vf_(&r2[nb + 1]) / vg_(&r2[nb + 1])) * (vf_(r2) /
		     vg_(r2));
	    svv += exp(st * 2);
	    if (l1 % 100 == 0 && twrite) {
		s_wsle(&io___68);
		r__1 = dt * l1;
		do_lio(&c__4, &c__1, (char *)&r__1, (ftnlen)sizeof(real));
		r__2 = sv / l1;
		do_lio(&c__4, &c__1, (char *)&r__2, (ftnlen)sizeof(real));
		e_wsle();
	    }

/* Konfiguration auf file schreiben */

	    if (cwrite) {
		i__3 = nb + 1;
		for (i = 0; i <= i__3; ++i) {
		    s_wsfe(&io___69);
		    for (j = 1; j <= 3; ++j) {
			do_fio(&c__1, (char *)&x[j + i * 3 - 1], (ftnlen)
				sizeof(real));
		    }
		    e_wsfe();
/* L300: */
		}
	    }
/* L200: */
	}

/* hier endet die Simulationsschleife ueber nc Konfigurationen */

/* berechne mittleren potentialanteil der Wirkung */
/* fuer eine anfangskonfiguration */
	sv /= nc;
/* Computing 2nd power */
	r__1 = sv;
	svv = svv / nc - r__1 * r__1;
	sa += sv;
	sav += svv;
/* L1000: */
    }

/* hier endet die Simulationschleife ueber na Anfangs */
/* konfigurationen */
    sa /= na;
    sav /= na;
    s_wsle(&io___70);
    r__1 = (nb + 1) * par_1.db;
    do_lio(&c__4, &c__1, (char *)&r__1, (ftnlen)sizeof(real));
    r__2 = 1 / ((nb + 1) * par_1.db);
    do_lio(&c__4, &c__1, (char *)&r__2, (ftnlen)sizeof(real));
    do_lio(&c__4, &c__1, (char *)&sa, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___71);
    r__1 = (nb + 1) * par_1.db;
    do_lio(&c__4, &c__1, (char *)&r__1, (ftnlen)sizeof(real));
    r__2 = 1 / ((nb + 1) * par_1.db);
    do_lio(&c__4, &c__1, (char *)&r__2, (ftnlen)sizeof(real));
    do_lio(&c__4, &c__1, (char *)&sav, (ftnlen)sizeof(real));
    r__3 = sqrt(sav);
    do_lio(&c__4, &c__1, (char *)&r__3, (ftnlen)sizeof(real));
    e_wsle();
    goto L1;
L9999:
    return 0;
} /* MAIN__ */

/* Subroutine */ int s_(integer *i, real *x, real *q, real *z, real *zq, 
	integer *nb, real *r2, real *r2t)
{
    /* System generated locals */
    real r__1, r__2;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer j;
    static real r, vp, rt;


/* berechnet die ableitung der Wirkung am iten Schritt */

    /* Parameter adjustments */
    --zq;
    --z;
    --q;
    --x;

    /* Function Body */
    r = 0.f;
    rt = 0.f;
    for (j = 1; j <= 3; ++j) {
/* Computing 2nd power */
	r__1 = x[j + *i * 3] - q[j + *i * 3];
	r += r__1 * r__1;
/* Computing 2nd power */
	r__1 = x[j + *i * 3];
	rt += r__1 * r__1;
/* L10: */
    }
    r2[*i] = r;
    r2t[*i] = rt;
    vp = -2.f / par_1.w * par_1.v0 * exp(-r / par_1.w) * par_1.db;
    if (*i == 0) {
	for (j = 1; j <= 3; ++j) {
/* Projectile: kinetic */
/* Computing 2nd power */
	    r__1 = par_1.b;
	    z[j] = par_1.mb * (x[j + (*i + 1) * 3] - x[j + *i * 3]) - vp * (x[
		    j + *i * 3] - q[j + *i * 3]) / 2.f - (x[j + *i * 3] - q[j 
		    + *i * 3]) / (r__1 * r__1);
/*             potential */
/*             endpoint */
/* Target: kinetic */
/* Computing 2nd power */
	    r__1 = par_1.xo;
/* Computing 2nd power */
	    r__2 = par_1.b;
	    zq[j] = par_1.f1 * (q[j + (*i + 1) * 3] - par_1.f2 * q[j + *i * 3]
		    ) - vp * (q[j + *i * 3] - x[j + *i * 3]) / 2.f - r__1 * 
		    r__1 * q[j + *i * 3] - (q[j + *i * 3] - x[j + *i * 3]) / (
		    r__2 * r__2);
/*         potential */
/*         endpoint */
/* L100: */
	}
    } else if (*i == *nb + 1) {
	for (j = 1; j <= 3; ++j) {
/* Projectile: kinetic */
/* Computing 2nd power */
	    r__1 = par_1.b;
	    z[j] = par_1.mb * (-x[j + *i * 3] + x[j + (*i - 1) * 3]) - vp * (
		    x[j + *i * 3] - q[j + *i * 3]) / 2.f - (x[j + *i * 3] - q[
		    j + *i * 3]) / (r__1 * r__1);
/*             potential */
/*             endpoint */
/* Target: kinetic */
/* Computing 2nd power */
	    r__1 = par_1.xo;
/* Computing 2nd power */
	    r__2 = par_1.b;
	    zq[j] = par_1.f1 * (-par_1.f2 * q[j + *i * 3] + q[j + (*i - 1) * 
		    3]) - vp * (q[j + *i * 3] - x[j + *i * 3]) / 2.f - r__1 * 
		    r__1 * q[j + *i * 3] - (q[j + *i * 3] - x[j + *i * 3]) / (
		    r__2 * r__2);
/*         potential */
/*         endpoint */
/* L110: */
	}
    } else {
	for (j = 1; j <= 3; ++j) {
/* Projectile: kinetic */
	    z[j] = par_1.mb * (x[j + (*i + 1) * 3] - x[j + *i * 3] * 2.f + x[
		    j + (*i - 1) * 3]) - vp * (x[j + *i * 3] - q[j + *i * 3]);
/*             potential */
/* Target: kinetic */
	    zq[j] = par_1.f1 * (q[j + (*i + 1) * 3] - par_1.f2 * 2.f * q[j + *
		    i * 3] + q[j + (*i - 1) * 3]) - vp * (q[j + *i * 3] - x[j 
		    + *i * 3]);
/*         potential */
/* L120: */
	}
    }
    return 0;
} /* s_ */

real v_(real *r)
{
    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double exp(doublereal);


/* potential (radialsymetrisch) */
/* r ist der radius quadrat !!! */

    ret_val = par_1.v0 * exp(-(*r) / par_1.w);
    return ret_val;
} /* v_ */

real vf_(real *r)
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static real alpha;


/* potential (radialsymetrisch) gefaltet mit grundzustand */
/* des harmonischen oszi */
/* r ist der radius quadrat !!! */

    alpha = par_1.w * par_1.mo * par_1.omega;
/* Computing 3rd power */
    r__1 = sqrt(alpha / (alpha + 1)), r__2 = r__1;
    ret_val = par_1.v0 * (r__2 * (r__1 * r__1)) * exp(-(*r) / par_1.w * alpha 
	    / (alpha + 1) * .97f);
    return ret_val;
} /* vf_ */

real vg_(real *r)
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static real alpha;


/* potential (radialsymetrisch) gefaltet mit grundzustand */
/* des harmonischen oszi */
/* r ist der radius quadrat !!! */

    alpha = par_1.w * par_1.mo * par_1.omega;
/* Computing 3rd power */
    r__1 = sqrt(alpha / (alpha + 1)), r__2 = r__1;
    ret_val = par_1.v0 * (r__2 * (r__1 * r__1)) * exp(-(*r) / par_1.w * alpha 
	    / (alpha + 1));
    return ret_val;
} /* vg_ */

real ran_(void)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1;

    /* Builtin functions */
    double d_int(doublereal *);

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 */
/*                                                                      C 
*/
/*  THIS ROUTINE CALCULATES A RANDOM NUMBER.                            C 
*/
/*                                                                      C 
*/
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 */
    rg_2.xk += .5;
    d__1 = rg_2.xk * 13136923.;
    rg_2.x = rg_2.x * 14662125. + d_int(&d__1);
    d__1 = rg_2.x / 281474976710656.;
    rg_2.x -= d_int(&d__1) * 281474976710656.;
    ret_val = rg_2.x / 281474976710656.;
    return ret_val;
} /* ran_ */

real rgau_(void)
{
    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), cos(doublereal);

    /* Local variables */
    static real x1, x2;
    extern real ran_(void);

    x1 = ran_();
    x2 = ran_();
    ret_val = sqrt(log(x1) * -2.f) * cos(const_1.pi * 2.f * x2);
    return ret_val;
} /* rgau_ */

/* Subroutine */ int makran_(real *w, real *wq, integer *n, real *ds)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal), log(doublereal), sqrt(doublereal)
	    ;

    /* Local variables */
    static integer i, j;
    static real x1, x2, lg, cs, sn;
    extern real ran_(void);


/* maximale zahl der beta schritte */


/* felder fuer zufallszahlen */


/* macht die zufallszahlen fuer einen schritt */
/* fuer das projektil */

    /* Parameter adjustments */
    --wq;
    --w;

    /* Function Body */
    i__1 = *n;
    for (i = 0; i <= i__1; i += 2) {
	for (j = 1; j <= 3; ++j) {
L101:
	    x1 = ran_();
	    if (x1 < 1e-16f) {
		goto L101;
	    }
	    x2 = ran_();
	    cs = cos(const_1.pi * 2.f * x2);
	    sn = sin(const_1.pi * 2.f * x2);
	    lg = sqrt(log(x1) * -2.f) * *ds;
	    w[j + i * 3] = lg * cs;
	    w[j + (i + 1) * 3] = lg * sn;
/* L100: */
	}
    }

/* zufallszahlen fuer harmonischen oszillator */

    i__1 = *n;
    for (i = 0; i <= i__1; i += 2) {
	for (j = 1; j <= 3; ++j) {
L111:
	    x1 = ran_();
	    if (x1 < 1e-16f) {
		goto L111;
	    }
	    x2 = ran_();
	    cs = cos(const_1.pi * 2.f * x2);
	    sn = sin(const_1.pi * 2.f * x2);
	    lg = sqrt(log(x1) * -2.f) * *ds;
	    wq[j + i * 3] = lg * cs;
	    wq[j + (i + 1) * 3] = lg * sn;
/* L110: */
	}
    }
    return 0;
} /* makran_ */

/* Main program alias */ int langv2_ () { MAIN__ (); return 0; }
