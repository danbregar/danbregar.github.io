/*
                      FoilSim II  - Airfoil  mode
                          veloicty effects only
   
                           A Java Applet
               to perform Kutta-Joukowski Airfoil analysis

                     Version 1.3b   - 16 Oct 00

                         Written by Tom Benson
                       NASA Glenn Research Center

>                              NOTICE
>This software is in the Public Domain.  It may be freely copied and used in
>non-commercial products, assuming proper credit to the author is given.  IT
>MAY NOT BE RESOLD.  If you want to use the software for commercial
>products, contact the author.
>No copyright is claimed in the United States under Title 17, U. S. Code.
>This software is provided "as is" without any warranty of any kind, either
>express, implied, or statutory, including, but not limited to, any warranty
>that the software will conform to specifications, any implied warranties of
>merchantability, fitness for a particular purpose, and freedom from
>infringement, and any warranty that the documentation will conform to the
>program, or any warranty that the software will be error free.
>In no event shall NASA be liable for any damages, including, but not
>limited to direct, indirect, special or consequential damages, arising out
>of, resulting from, or in any way connected with this software, whether or
>not based on warranty, contract, tort or otherwise, whether or not injury
>was sustained by persons or property or otherwise, and whether or not loss
>was sustained from, or arose out of the results of, or use of, the software
>or services provided hereunder.
 
  New test -
             keep most of FoilSim but only use the velocity input
 
                                           TJB  16 Oct 00

*/

import java.awt.*;
import java.lang.Math ;

public class Foilvel extends java.applet.Applet {
 
   static double convdr = 3.1415926/180. ;
   static double pid2 = 3.1415926/2.0 ;
   static double rval,ycval,xcval,gamval,alfval,thkval,camval,chrd,cl ;
   static double thkinpt,caminpt ;                 /* MODS 10 Sep 99 */
   static double leg,teg,lem,tem;
   static double usq,vsq,alt,altmax,area,armax,armin ;
   static double chord,span,aspr,arold,chrdold,spnold ; /* Mod 13 Jan 00 */
   static double q0,ps0,pt0,ts0,rho,rlhum,temf,presm ;
   static double lyg,lrg,lthg,lxgt,lygt,lrgt,lthgt;/* MOD 20 Jul */
   static double lxm,lym,lxmt,lymt,vxdir;/* MOD 20 Jul */
   static double deltb,xflow ;             /* MODS  20 Jul 99 */
   static double delx,delt,vfsd,spin,spindr,yoff ;
   static double vel,pres,lift,side,omega,radcrv,relsy,angr ;

   static double rg[][]  = new double[20][40] ; 
   static double thg[][] = new double[20][40] ; 
   static double xg[][]  = new double[20][40] ; 
   static double yg[][]  = new double[20][40] ; 
   static double xm[][]  = new double[20][40] ; 
   static double ym[][]  = new double[20][40] ; 
   static double plp[]   = new double[40] ;
   static double plv[]   = new double[40] ;

   int nptc,npt2,nlnc,nln2,rdflag,browflag,probflag,anflag;
   int foil,flflag,lunits,lftout,planet ;
   int conflag,displ,dispp,antim,ancol;   /* MODS  2 APR 99 - 22 APR -27 JUL */
       /* units data */
   static double vmn,almn,angmn,vmx,almx,angmx ;
   static double camn,thkmn,camx,thkmx ;
   static double chrdmn,spnmn,armn,chrdmx,spnmx,armx ;
   static double vconv,vmaxa,vmaxb ;
   static double pconv,pmax,pmin,lconv,fconv,fmax,fmaxb;
   int lflag,gflag,plscale,nond;
       /*  plot & probe data */
   static double fact,xpval,ypval,pbval,factp;
   static double prg,pthg,pxg,pyg,pxm,pym ;
   int pboflag,xt,yt,ntikx,ntiky,npt,xtp,ytp ;
   int lines,nord,nabs,ntr ;
   static double begx,endx,begy,endy ;
   static String labx,labxu,laby,labyu ;
   static double pltx[][]  = new double[3][40] ;
   static double plty[][]  = new double[3][40] ;

   Solver solve ;
   Slvpnl slvpnl ;
   L l ;
   CardLayout layin ;
   Image offImg1 ;
   Graphics off1Gg ;
   Image offImg2 ;
   Graphics off2Gg ;
   Image offImg3 ;
   Graphics off3Gg ;

   public void init() {
     int i;
     Foilvel a = new Foilvel() ;
     solve = new Solver() ;

     offImg1 = createImage(this.size().width,
                      this.size().height) ;
     off1Gg = offImg1.getGraphics() ;
     offImg2 = createImage(this.size().width,
                      this.size().height) ;
     off2Gg = offImg2.getGraphics() ;
     offImg3 = createImage(this.size().width,
                      this.size().height) ;
     off3Gg = offImg3.getGraphics() ;

     setLayout(new GridLayout(1,2,5,5)) ;

     solve.setDefaults () ;
 
     slvpnl = new Slvpnl(this) ;
     l = new L(this) ;

     add(l) ;
     add(slvpnl) ;

     solve.getFreeStream ();
     computeFlow () ;
     l.view.start() ;
     slvpnl.outpnl.start() ;
  }
 
  public Insets insets() {
     return new Insets(5,5,5,5) ;
  }

  public void computeFlow() { 

     if (flflag == 1) {
         solve.getCirc ();                   /* get circulation */
         solve.genFlow () ;
         solve.getFreeStream () ;
     }
 
     loadOut() ;
     slvpnl.outpnl.loadPlot() ;
  }

  public int filter0(double inumbr) {
        //  output only to .
       int number ;
       int intermed ;
 
       number = (int) (inumbr);
       return number ;
  }

  public float filter3(double inumbr) {
     //  output only to .001
       float number ;
       int intermed ;
 
       intermed = (int) (inumbr * 1000.) ;
       number = (float) (intermed / 1000. );
       return number ;
  }
 
  public float filter5(double inumbr) {
     //  output only to .00001
       float number ;
       int intermed ;
 
       intermed = (int) (inumbr * 100000.) ;
       number = (float) (intermed / 100000. );
       return number ;
  }
 
  public void setUnits() {   // Switching Units
       double ovs,chords,spans,aros,chos,spos ;
       double alts,ares ;

       alts = alt / lconv ;
       chords = chord / lconv ;
       spans = span / lconv ;
       ares = area /lconv/lconv ;
       aros = arold /lconv/lconv ;
       chos = chrdold / lconv ;
       spos = spnold / lconv ;
       ovs = vfsd / vconv ;

       switch (lunits) {
          case 0: {                             /* English */
            lconv = 1.;                      /*  feet    */
            vconv = .6818; vmaxa = 250.; vmaxb = 100. ;  /*  mph  */
            fconv = 1.0; fmax = 100000.; fmaxb = .5;  /* pounds   */
            pconv = 14.7  ;                   /* lb/sq in */
            break;
          }
          case 1: {                             /* Metric */
            lconv = .3048;                    /* meters */
            vconv = 1.097;  vmaxa = 400.; vmaxb = 167.;   /* km/hr  */
            fconv = 4.448 ; fmax = 500000.; fmaxb = 2.5; /* newtons */
            pconv = 101.3 ;               /* kilo-pascals */
            break ;
          }
       }
 
       alt = alts * lconv ;
       chord = chords * lconv ;
       span = spans * lconv ;
       area = ares * lconv * lconv ;
       arold = aros * lconv * lconv ;
       chrdold = chos * lconv ;
       spnold = spos * lconv ;
       vfsd  = ovs * vconv;

       return ;
  }

  public void loadOut() {   // output routine
     double stfact ;
                          // stall model
     stfact = 1.0 ;
     if (anflag == 1) {
         if (alfval > 10.0 ) {
            stfact = .5 + .1 * alfval - .005 * alfval * alfval ;
         }
         if (alfval < -10.0 ) {
            stfact = .5 - .1 * alfval - .005 * alfval * alfval ;
         }
     }

     if (lftout == 1) {
       lift = cl * stfact ;
       slvpnl.conpnl.out.setText(String.valueOf(filter3(lift))) ;
     }
     if (lftout == 0) {
       lift = cl * stfact * q0 * area / lconv / lconv ;    /* lift in lbs */
       lift = lift * fconv ;
       if (Math.abs(lift) <= 10.0) {
          slvpnl.conpnl.out.setText(String.valueOf(filter3(lift))) ;
       }
       if (Math.abs(lift) > 10.0) {
          slvpnl.conpnl.out.setText(String.valueOf(filter0(lift))) ;
       }
     }
 
     return ;
  }

  class Solver {
 
     Solver () {
     }

     public void setDefaults() {

        planet = 0 ;
        lunits = 0 ;
        lftout = 0 ;
        nlnc = 15 ;
        nln2 = nlnc/2 + 1 ;
        nptc = 37 ;
        npt2 = nptc/2 + 1 ;
        deltb = .5 ;
        foil = 1 ;
        flflag = 1;
        thkval = .5 ;
        thkinpt = 12.5 ;                   /* MODS 10 SEP 99 */
        camval = 0.0 ;
        caminpt = 0.0 ;
        alfval = 5.0 ;
        gamval = 0.0 ;
        spin = 0.0 ;
        spindr = 1.0 ;
        rval = 1.0 ;
        ycval = 0.0 ;
        xcval = 0.0 ;
        conflag = 1 ;                             /* MODS  2 Apr 99 */
        displ = 0 ;                              /* MODS  22 Apr 99 */
        dispp = 5 ;
        lift = 313. ;
 
        xpval = 2.1;
        ypval = -.5 ;
        pboflag = 0 ;
        xflow = -10.0;                             /* MODS  20 Jul 99 */

        pconv = 14.7;
        pmin = .5 ;
        pmax = 1.0 ;
        fconv = 1.0 ;
        fmax = 100000. ;
        fmaxb = .50 ;
        vconv = .6818 ;
        vfsd = 100. ;
        vmaxa = 250. ;
        vmaxb = 100. ;
        lconv = 1.0 ;

        alt = 0.0 ;
        altmax = 50000. ;
        chrdold = chord = 2.0 ;
        spnold = span = 10.0 ;
        aspr = 5.0 ;
        arold = area = 20.0 ;
        armax = 1000.01 ;
        armin = .01 ;                 /* MODS 9 SEP 99 */
 
        xt = 120;  yt = 75; fact = 30.0 ;
        xtp = 95; ytp = 130; factp = 25.0 ;
 
        probflag = 2 ;
        anflag = 0 ;
        vmn = 0.0;     vmx = 250.0 ;
        almn = 0.0;    almx = 50000.0 ;
        angmn = -20.0; angmx = 20.0 ;
        camn = -25.0;  camx = 25.0 ;
        thkmn = 1.0; thkmx = 26.0 ;
        chrdmn = .1 ;  chrdmx = 10.1 ;
        spnmn = .1 ;  spnmx = 100.1 ;
        armn = .01 ;  armx = 1000.01 ;

        return ;
     }

     public void getFreeStream() {    //  free stream conditions
       double hite,pvap,rgas,gama ;       /* MODS  19 Jan 00  whole routine*/

       rgas = 1718. ;                /* ft2/sec2 R */
       gama = 1.4 ;
       hite = alt/lconv ;
       if (planet == 0) {    // Earth  standard day
         if (conflag == 1) {
           if (hite <= 36152.) {           // Troposphere
              ts0 = 518.6 - 3.56 * hite/1000. ;
              ps0 = 2116. * Math.pow(ts0/518.6,5.256) ;
           }
           if (hite >= 36152. && hite <= 82345.) {   // Stratosphere
              ts0 = 389.98 ;
              ps0 = 2116. * .2236 *
                 Math.exp((36000.-hite)/(53.35*389.98)) ;
           }
           if (hite >= 82345.) {
              ts0 = 389.98 + 1.645 * (hite-82345)/1000. ;
              ps0 = 2116. *.02456 * Math.pow(ts0/389.98,-11.388) ;
           }
           rlhum = 0.0 ;
           temf = ts0 - 459.6 ;
           if (temf <= 0.0) temf = 0.0 ;                    
           presm = ps0 * 29.92 / 2116. ;
         }
         if (conflag == 2) {
            ts0 = temf + 459.6 ;
            if (temf < 0.0) {
                  temf = 0.0 ;
                  rlhum = 0.0 ;
            }
             ps0 = presm * 2116. / 29.92 ;
         }
         pvap = rlhum*(2.685+.00353*Math.pow(temf,2.245));/* Eq 1:6A  Domasch */
         rho = (ps0 - .379*pvap)/(rgas * ts0) ;  /* effect of humidty */
         rho = ps0/(53.3 * 32.17 * ts0) ;
       }

       if (planet == 1) {   // Mars - curve fit of orbiter data
         rgas = 1149. ;                /* ft2/sec2 R */
         gama = 1.29 ;

         if (hite <= 22960.) {
            ts0 = 434.02 - .548 * hite/1000. ;
            ps0 = 14.62 * Math.pow(2.71828,-.00003 * hite) ;
         }
         if (hite > 22960.) {
            ts0 = 449.36 - 1.217 * hite/1000. ;
            ps0 = 14.62 * Math.pow(2.71828,-.00003 * hite) ;
         }
         rho = ps0/(rgas*ts0) ;
       }

       q0  = .5 * rho * vfsd * vfsd / (vconv * vconv) ;
       pt0 = ps0 + q0 ;

       return ;
     }

     public void getCirc() {   // circulation from Kutta condition
       double thet,rdm,thtm ;
       double beta,rball;
       int index;

       xcval = 0.0 ;
       switch (foil)  {
          case 0: {         /* get circulation from spin for baseball */
              rball = .1 ;         /* baseball radius = .1 ft = 1.2 in */
              gamval = 4.0 * 3.1415926 * 3.1415926 *spin * rball * rball
                                 / (vfsd/vconv) ;
              gamval = gamval * spindr ;
              ycval = .0001 ;
              break ;
          }
          case 1:  {                  /* Juokowski geometry*/
              ycval = camval / 2.0 ;
              rval = thkval/4.0 +Math.sqrt(thkval*thkval/16.0+ycval*ycval +1.0);
              xcval = 1.0 - Math.sqrt(rval*rval - ycval*ycval) ;
              beta = Math.asin(ycval/rval)/convdr ;     /* Kutta condition */
              gamval = 2.0*rval*Math.sin((alfval+beta)*convdr) ;
              break ;
          }
          case 2:  {                  /* Elliptical geometry*/
              ycval = camval / 2.0 ;
              rval = thkval/4.0 + Math.sqrt(thkval*thkval/16.0+ycval*ycval+1.0);
              beta = Math.asin(ycval/rval)/convdr ;    /* Kutta condition */
              gamval = 2.0*rval*Math.sin((alfval+beta)*convdr) ;
              break ;
          }
          case 3:  {                  /* Plate geometry*/
              ycval = camval / 2.0 ;
              rval = Math.sqrt(ycval*ycval+1.0);
              beta = Math.asin(ycval/rval)/convdr ;    /* Kutta condition */
              gamval = 2.0*rval*Math.sin((alfval+beta)*convdr) ;
              break ;
          }
       }
                                                   /* geometry */
       for (index =1; index <= nptc; ++index) {
           thet = (index -1)*360./(nptc-1) ;
           xg[0][index] = rval * Math.cos(convdr * thet) + xcval ;
           yg[0][index] = rval * Math.sin(convdr * thet) + ycval ;
           rg[0][index] = Math.sqrt(xg[0][index]*xg[0][index] +
                                    yg[0][index]*yg[0][index])  ;
           thg[0][index] = Math.atan2(yg[0][index],xg[0][index])/convdr;
           xm[0][index] = (rg[0][index] + 1.0/rg[0][index])*
                        Math.cos(convdr*thg[0][index]) ;
           ym[0][index] = (rg[0][index] - 1.0/rg[0][index])*
                        Math.sin(convdr*thg[0][index]) ;
           rdm = Math.sqrt(xm[0][index]*xm[0][index] +
                           ym[0][index]*ym[0][index])  ;
           thtm = Math.atan2(ym[0][index],xm[0][index])/convdr;
           xm[0][index] = rdm * Math.cos((thtm - alfval)*convdr);
           ym[0][index] = rdm * Math.sin((thtm - alfval)*convdr);
           getVel(rval,thet) ;
           plp[index] = ((ps0 + pres * q0)/2116.) * pconv ;
           plv[index] = vel * vfsd ;
       }

       return ;
     }

     public void genFlow() {   // generate flowfield
       double rnew,thet,psv,fxg;
       int k,index;
                              /* all lines of flow  except stagnation line*/
       for (k=1; k<=nlnc; ++k) {
         psv = -.5*(nln2-1) + .5*(k-1) ;
         fxg = xflow ;
         for (index =1; index <=nptc; ++ index) {
           solve.getPoints (fxg,psv) ;
           xg[k][index]  = lxgt ;
           yg[k][index]  = lygt ;
           rg[k][index]  = lrgt ;
           thg[k][index] = lthgt ;
           xm[k][index]  = lxmt ;
           ym[k][index]  = lymt ;
           if (anflag == 1) {           // stall model
              if (alfval > 10.0 && psv > 0.0) {
                   if (xm[k][index] > 0.0) {
                      ym[k][index] = ym[k][index -1] ;
                   }
              }
              if (alfval < -10.0 && psv < 0.0) {
                   if (xm[k][index] > 0.0) {
                      ym[k][index] = ym[k][index -1] ;
                   }
              }
           }
           solve.getVel(lrg,lthg) ;
           fxg = fxg + vxdir*deltb ;
         }
       }
                                              /*  stagnation line */
       k = nln2 ;
       psv = 0.0 ;
                                              /*  incoming flow */
       for (index =1; index <= npt2; ++ index) {
           rnew = 10.0 - (10.0 - rval)*Math.sin(pid2*(index-1)/(npt2-1)) ;
           thet = Math.asin(.999*(psv - gamval*Math.log(rnew/rval))/
                                   (rnew - rval*rval/rnew)) ;
           fxg =  - rnew * Math.cos(thet) ;
           solve.getPoints (fxg,psv) ;
           xg[k][index]  = lxgt ;
           yg[k][index]  = lygt ;
           rg[k][index]  = lrgt ;
           thg[k][index] = lthgt ;
           xm[k][index]  = lxmt ;
           ym[k][index]  = lymt ;
       }
                                              /*  downstream flow */
       for (index = 1; index <= npt2; ++ index) {
           rnew = 10.0 + .01 - (10.0 - rval)*Math.cos(pid2*(index-1)/(npt2-1)) ;
           thet = Math.asin(.999*(psv - gamval*Math.log(rnew/rval))/
                                      (rnew - rval*rval/rnew)) ;
           fxg =   rnew * Math.cos(thet) ;
           solve.getPoints (fxg,psv) ;
           xg[k][npt2+index]  = lxgt ;
           yg[k][npt2+index]  = lygt ;
           rg[k][npt2+index]  = lrgt ;
           thg[k][npt2+index] = lthgt ;
           xm[k][npt2+index]  = lxmt ;
           ym[k][npt2+index]  = lymt ;
       }
                                              /*  stagnation point */
       xg[k][npt2]  = xcval ;
       yg[k][npt2]  = ycval ;
       rg[k][npt2]  = Math.sqrt(xcval*xcval+ycval*ycval) ;
       thg[k][npt2] = Math.atan2(ycval,xcval)/convdr ;
       xm[k][npt2]  = (xm[k][npt2+1] + xm[k][npt2-1])/2.0 ;
       ym[k][npt2]  = (ym[0][nptc/4+1] + ym[0][nptc/4*3+1])/2.0 ;
                                /*  compute lift coefficient */
       leg = xcval - Math.sqrt(rval*rval - ycval*ycval) ;
       teg = xcval + Math.sqrt(rval*rval - ycval*ycval) ;
       lem = leg + 1.0/leg ;
       tem = teg + 1.0/teg ;
       chrd = tem - lem ;
       cl = gamval*4.0*3.1415926/chrd ;

       return ;
     }

     public void getPoints(double fxg, double psv) {   // flow in x-psi
       double radm,thetm ;                /* MODS  20 Jul 99  whole routine*/
       double fnew,ynew,yold,rfac,deriv ;
       double xold,xnew,thet ;
       double rmin,rmax ;
       int iter,isign;
                       /* get variables in the generating plane */
                           /* iterate to find value of yg */
       ynew = 10.0 ;
       yold = 10.0 ;
       if (psv < 0.0) ynew = -10.0 ;
       if (Math.abs(psv) < .001 && alfval < 0.0) ynew = rval ;
       if (Math.abs(psv) < .001 && alfval >= 0.0) ynew = -rval ;
       fnew = 0.1 ;
       iter = 1 ;
       while (Math.abs(fnew) >= .00001 && iter < 25) {
           ++iter ;
           rfac = fxg*fxg + ynew*ynew ;
           if (rfac < rval*rval) rfac = rval*rval + .01 ;
           fnew = psv - ynew*(1.0 - rval*rval/rfac)
                  - gamval*Math.log(Math.sqrt(rfac)/rval) ;
           deriv = - (1.0 - rval*rval/rfac)
               - 2.0 * ynew*ynew*rval*rval/(rfac*rfac)
               - gamval * ynew / rfac ;
           yold = ynew ;
           ynew = yold  - .5*fnew/deriv ;
       }
       lyg = yold ;
                                     /* rotate for angle of attack */
       lrg = Math.sqrt(fxg*fxg + lyg*lyg) ;
       lthg = Math.atan2(lyg,fxg)/convdr ;
       lxgt = lrg * Math.cos(convdr*(lthg + alfval)) ;
       lygt = lrg * Math.sin(convdr*(lthg + alfval)) ;
                              /* translate cylinder to generate airfoil */
       lxgt = lxgt + xcval ;
       lygt = lygt + ycval ;
       lrgt = Math.sqrt(lxgt*lxgt + lygt*lygt) ;
       lthgt = Math.atan2(lygt,lxgt)/convdr ;
                               /*  Kutta-Joukowski mapping */
       lxm = (lrgt + 1.0/lrgt)*Math.cos(convdr*lthgt) ;
       lym = (lrgt - 1.0/lrgt)*Math.sin(convdr*lthgt) ;
                              /* tranforms for view fixed with free stream */
                /* take out rotation for angle of attack mapped and cylinder */
       radm = Math.sqrt(lxm*lxm+lym*lym) ;
       thetm = Math.atan2(lym,lxm)/convdr ;
       lxmt = radm*Math.cos(convdr*(thetm-alfval)) ;
       lymt = radm*Math.sin(convdr*(thetm-alfval)) ;

       lxgt = lxgt - xcval ;
       lygt = lygt - ycval ;
       lrgt = Math.sqrt(lxgt*lxgt + lygt*lygt)  ;
       lthgt = Math.atan2(lygt,lxgt)/convdr;
       lxgt = lrgt * Math.cos((lthgt - alfval)*convdr);
       lygt = lrgt * Math.sin((lthgt - alfval)*convdr);

       return ;
     }
 
     public void getVel(double radius, double theta) {  //velocity and pressure 
      double ur,uth,jake1,jake2,jakesq ;
      double xloc,yloc,thrad,alfrad ;

      thrad = convdr * theta ;
      alfrad = convdr * alfval ;
                                /* get x, y location in cylinder plane */
      xloc = radius * Math.cos(thrad) ;
      yloc = radius * Math.sin(thrad) ;
                                /* velocity in cylinder plane */
      ur  = Math.cos(thrad-alfrad)*(1.0-(rval*rval)/(radius*radius)) ;
      uth = -Math.sin(thrad-alfrad)*(1.0+(rval*rval)/(radius*radius))
                            - gamval/radius;
      usq = ur*ur + uth*uth ;
      vxdir = ur * Math.cos(thrad) - uth * Math.sin(thrad) ; // MODS  20 Jul 99 
                                /* translate to generate airfoil  */
      xloc = xloc + xcval ;
      yloc = yloc + ycval ;
                                   /* compute new radius-theta  */
      radius = Math.sqrt(xloc*xloc + yloc*yloc) ;
      thrad  = Math.atan2(yloc,xloc) ;
                                   /* compute Joukowski Jacobian  */
      jake1 = 1.0 - Math.cos(2.0*thrad)/(radius*radius) ;
      jake2 = Math.sin(2.0*thrad)/(radius*radius) ;
      jakesq = jake1*jake1 + jake2*jake2 ;
      if (Math.abs(jakesq) <= .01) jakesq = .01 ;  /* protection */
      vsq = usq / jakesq ;
          /* vel is velocity ratio - pres is coefficient  (p-p0)/q0   */
      if (foil > 0) {
           vel = Math.sqrt(vsq) ;
           pres = 1.0 - vsq ;
      }
      if (foil == 0) {
           vel = Math.sqrt(usq) ;
           pres = 1.0 - usq ;
      }
      return ;
    }

  }
 
  class L extends Panel {
     Foilvel outerparent ;
     Viewer view ;
     Inppnl inppnl ;

     L (Foilvel target) { 
         outerparent = target ;
         setLayout(new GridLayout(2,1,5,5)) ;

         view  = new Viewer(outerparent) ;
         inppnl  = new Inppnl(outerparent) ;

         add(view) ;
         add(inppnl) ;
     }

     class Inppnl extends Panel {
        Foilvel outerparent ;
        Lft lft ;
        Rght rght ;

        Inppnl (Foilvel target) {

           outerparent = target ;
           setLayout(new GridLayout(1,2,5,5)) ;

           lft = new Lft(outerparent) ;
           rght = new Rght(outerparent) ;

           add(lft) ;
           add(rght) ;
        }

        class Lft extends Panel {
           Foilvel outerparent ;
           TextField f1,f2,o1,o2 ;
           Label l1,l2 ;
           Label l01,l02 ;
           Label lo1,lo2 ;
      
           Lft (Foilvel target) {
      
             outerparent = target ;
             setLayout(new GridLayout(7,2,2,10)) ;

             l1 = new Label("Speed-mph", Label.CENTER) ;
             f1 = new TextField("100.0",5) ;

             add(l1) ;
             add(f1) ;

             add(new Label("Change", Label.RIGHT)) ;
             add(new Label("value:", Label.LEFT)) ;

             add(new Label("1. Back", Label.RIGHT)) ;
             add(new Label("space over", Label.LEFT)) ;

             add(new Label("old", Label.RIGHT)) ;
             add(new Label("value.", Label.LEFT)) ;

             add(new Label("2. Enter", Label.RIGHT)) ;
             add(new Label("new value.", Label.LEFT)) ;

             add(new Label("3. Hit", Label.RIGHT)) ;
             add(new Label("ENTER key.", Label.LEFT)) ;

             add(new Label(" ", Label.RIGHT)) ;
             add(new Label(" ", Label.LEFT)) ;
           }

           public boolean handleEvent(Event evt) {
             Double V1,V2 ;
             double v1,v2 ;
             float fl1 ;
             int i1,i2 ;

             if(evt.id == Event.ACTION_EVENT) {
               V1 = Double.valueOf(f1.getText()) ;
               v1 = V1.doubleValue() ;

               vfsd = v1 ;
               if(v1 < vmn) {
                 vfsd = v1 = vmn ;
                 fl1 = (float) v1 ;
                 f1.setText(String.valueOf(fl1)) ;
               }
               if(v1 > vmx) {
                 vfsd = v1 = vmx ;
                 fl1 = (float) v1 ;
                 f1.setText(String.valueOf(fl1)) ;
               }

               i1 = (int) (((v1 - vmn)/(vmx-vmn))*1000.) ;
    
               rght.s1.setValue(i1) ;

               computeFlow() ;
               return true ;
             }
             else return false ;
           } // Handler
         }  // Inleft

         class Rght extends Panel {
            Foilvel outerparent ;
            Scrollbar s1,s2;
            Choice plntch;

            Rght (Foilvel target) {
             int i1,i2 ;

             outerparent = target ;
             setLayout(new GridLayout(7,1,2,10)) ;

             i1 = (int) (((100.0 - vmn)/(vmx-vmn))*1000.) ;

             s1 = new Scrollbar(Scrollbar.HORIZONTAL,i1,10,0,1000);

             add(s1) ;
             add(new Label(" ", Label.CENTER)) ;
             add(new Label("1. Move slider bar", Label.CENTER)) ;
             add(new Label("2. Click on arrows", Label.CENTER)) ;
             add(new Label(" ", Label.CENTER)) ;
             add(new Label(" ", Label.CENTER)) ;
             add(new Label(" ", Label.CENTER)) ;

           }

           public boolean handleEvent(Event evt) {
                if(evt.id == Event.ACTION_EVENT) {
                   this.handleBar(evt) ;
                   return true ;
                }
                if(evt.id == Event.SCROLL_ABSOLUTE) {
                   this.handleBar(evt) ;
                   return true ;
                }
                if(evt.id == Event.SCROLL_LINE_DOWN) {
                   this.handleBar(evt) ;
                   return true ;
                }
                if(evt.id == Event.SCROLL_LINE_UP) {
                   this.handleBar(evt) ;
                   return true ;
                }
                if(evt.id == Event.SCROLL_PAGE_DOWN) {
                   this.handleBar(evt) ;
                   return true ;
                }
                if(evt.id == Event.SCROLL_PAGE_UP) {
                   this.handleBar(evt) ;
                   return true ;
                }
                else return false ;
           }

           public void handleBar(Event evt) {
              int i1,i2 ;
              double v1,v2 ;
              float fl1,fl2 ;

     // Input for computations
              i1 = s1.getValue() ;

              vfsd   = v1 = i1 * (vmx - vmn)/ 1000. + vmn ;

              fl1 = (float) v1 ;

              lft.f1.setText(String.valueOf(fl1)) ;
       
              computeFlow() ;
            }
         }  // Inright 
     }  // Inppnl 

     class Viewer extends Canvas  
         implements Runnable{
        Foilvel outerparent ;
        Thread runner ;
        Point locate,anchor;
   
        Viewer (Foilvel target) {
            setBackground(Color.black) ;
            runner = null ;
        }    

        public Insets insets() {
           return new Insets(0,10,0,10) ;
        }
 
        public void start() {
           if (runner == null) {
              runner = new Thread(this) ;
              runner.start() ;
           }
           antim = 0 ;                              /* MODS  21 JUL 99 */
           ancol = 1 ;                              /* MODS  27 JUL 99 */
        }
   
        public void run() {
          int timer ;
    
          timer = 100 ;
          while (true) {
             ++ antim ;
             try { Thread.sleep(timer); }
             catch (InterruptedException e) {}
             l.view.repaint() ;
             if (antim == 3) {
                antim = 0;
                ancol = - ancol ;               /* MODS 27 JUL 99 */
             }
             timer = 135 - (int) (.227 *vfsd/vconv) ;
          }
        }
   
        public void update(Graphics g) {
           l.view.paint(g) ;
        }
    
        public void paint(Graphics g) {
           int i,j,k,n ;
           int xlabel,ylabel,ind,inmax,inmin ;
           int exes[] = new int[8] ;
           int whys[] = new int[8] ;
           double offx,scalex,offy,scaley,waste,incy,incx;
           double xtrans,ytrans,xl,yl;
           int camx[] = new int[19] ;
           int camy[] = new int[19] ;
           Color col ;
   
           col = new Color(0,0,0) ;
           if(planet == 0) col = Color.cyan ;
           if(planet == 1) col = Color.orange ;
           off1Gg.setColor(Color.black) ;
           off1Gg.fillRect(0,0,300,300) ;
   
           if (displ == 4) {              // Top View
             off1Gg.setColor(Color.white) ;
             exes[0] = (int) (.25*fact*(-span)) + xt ;
             whys[0] = (int) (.25*fact*(-chord)) + yt ;
             exes[1] = (int) (.25*fact*(-span)) + xt ;
             whys[1] = (int) (.25*fact*(chord)) + yt ;
             exes[2] = (int) (.25*fact*(span)) + xt ;
             whys[2] = (int) (.25*fact*(chord)) + yt ;
             exes[3] = (int) (.25*fact*(span)) + xt ;
             whys[3] = (int) (.25*fact*(-chord)) + yt ;
             off1Gg.fillPolygon(exes,whys,4) ;
             off1Gg.setColor(Color.green) ;
             off1Gg.drawString("View - Top",10,10) ;
             off1Gg.drawString("Span",exes[2]-20,whys[1]+20) ;
             off1Gg.drawLine(exes[0],whys[1]+5,exes[2],whys[1]+5) ;
             off1Gg.drawString("Chord",exes[2]+10,55) ;
             off1Gg.drawLine(exes[2]+5,whys[0],exes[2]+5,whys[1]) ;
             off1Gg.drawString("Flow",10,175) ;
             off1Gg.drawLine(20,165,20,135) ;
             exes[0] = 15 ;  exes[1] = 25; exes[2] = 20 ;
             whys[0] = 135 ;  whys[1] = 135; whys[2] = 125 ;
             off1Gg.fillPolygon(exes,whys,3) ;
           }

           if (displ <= 3) {  // Side View
            if (vfsd > .01) {
                                               /* plot airfoil flowfield */
             for (j=1; j<=nln2-1; ++j) {           /* lower half */
                exes[1] = (int) (fact*xm[j][1]) + xt ;
                whys[1] = (int) (fact*(-ym[j][1])) + yt ;
                for (i=2 ; i<= nptc; ++i) {
                   exes[0] = exes[1] ;
                   whys[0] = whys[1] ;
                   exes[1] = (int) (fact*xm[j][i]) + xt ;
                   whys[1] = (int) (fact*(-ym[j][i])) + yt ;
                   if (displ == 2) {                   /* MODS  21 JUL 99 */
                     off1Gg.setColor(Color.yellow) ;
                     off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                   }
                   if (displ == 1 && (i/3*3 == i) ) {
                     off1Gg.setColor(col) ;
                     for (n=1 ; n <= 4 ; ++n) {
                        if(i == 6 + (n-1)*9) off1Gg.setColor(Color.red);
                     }
                     if(i/9*9 == i) off1Gg.setColor(Color.white);
                     off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                   }
                   if (displ == 0 && ((i-antim)/3*3 == (i-antim)) ) {
                     if (ancol == -1) {          /* MODS  27 JUL 99 */
                       if((i-antim)/6*6 == (i-antim))off1Gg.setColor(col);
                       if((i-antim)/6*6 != (i-antim))off1Gg.setColor(Color.white);
                     }
                     if (ancol == 1) {          /* MODS  27 JUL 99 */
                       if((i-antim)/6*6 == (i-antim))off1Gg.setColor(Color.white);
                       if((i-antim)/6*6 != (i-antim))off1Gg.setColor(col);
                     }
                     off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                   }
                }
             }
             for (j=nln2+1; j<=nlnc; ++j) {          /* upper half */
                exes[1] = (int) (fact*xm[j][1]) + xt ;
                whys[1] = (int) (fact*(-ym[j][1])) + yt ;
                for (i=2 ; i<= nptc; ++i) {
                   exes[0] = exes[1] ;
                   whys[0] = whys[1] ;
                   exes[1] = (int) (fact*xm[j][i]) + xt ;
                   whys[1] = (int) (fact*(-ym[j][i])) + yt ;
                   if (displ == 2) {                     /* MODS  21 JUL 99 */
                     off1Gg.setColor(col) ;
                     off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                   }
                   if (displ == 1 && (i/3*3 == i) ) {
                     off1Gg.setColor(col);   /* MODS  27 JUL 99 */
                     for (n=1 ; n <= 4 ; ++n) {
                        if(i == 6 + (n-1)*9) off1Gg.setColor(Color.red);
                     }
                     if(i/9*9 == i) off1Gg.setColor(Color.white);
                     off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                   }
                   if (displ == 0 && ((i-antim)/3*3 == (i-antim)) ) {
                     if (ancol == -1) {          /* MODS  27 JUL 99 */
                       if((i-antim)/6*6 == (i-antim))off1Gg.setColor(col);
                       if((i-antim)/6*6 != (i-antim))off1Gg.setColor(Color.white);
                     }
                     if (ancol == 1) {          /* MODS  27 JUL 99 */
                       if((i-antim)/6*6 == (i-antim))off1Gg.setColor(Color.white);
                       if((i-antim)/6*6 != (i-antim))off1Gg.setColor(col);
                     }
                     off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                   }
                }
             }
             off1Gg.setColor(Color.white) ; /* stagnation */
             exes[1] = (int) (fact*xm[nln2][1]) + xt ;
             whys[1] = (int) (fact*(-ym[nln2][1])) + yt ;
             for (i=2 ; i<= npt2-1; ++i) {
                   exes[0] = exes[1] ;
                   whys[0] = whys[1] ;
                   exes[1] = (int) (fact*xm[nln2][i]) + xt ;
                   whys[1] = (int) (fact*(-ym[nln2][i])) + yt ;
                   if (displ <= 2) {             /* MODS  21 JUL 99 */
                     off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                   }
             }
             exes[1] = (int) (fact*xm[nln2][npt2+1]) + xt ;
             whys[1] = (int) (fact*(-ym[nln2][npt2+1])) + yt ;
             for (i=npt2+2 ; i<= nptc; ++i) {
                   exes[0] = exes[1] ;
                   whys[0] = whys[1] ;
                   exes[1] = (int) (fact*xm[nln2][i]) + xt ;
                   whys[1] = (int) (fact*(-ym[nln2][i])) + yt ;
                   if (displ <= 2) {                         /* MODS  21 JUL 99 */
                     off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                   }
             }
                                                  /*  probe location */
             if (pboflag > 0) {
                off1Gg.setColor(Color.red) ;
                off1Gg.fillOval((int) (fact*pxm) + xt,
                     (int) (fact*(-pym)) + yt - 2,5,5);
                off1Gg.setColor(Color.white) ;
                exes[0] = (int) (fact*(pxm + .1)) +xt ;
                whys[0] = (int) (fact*(-pym)) + yt ;
                exes[1] = (int) (fact*(pxm + .5)) +xt ;
                whys[1] = (int) (fact*(-pym)) + yt ;
                exes[2] = (int) (fact*(pxm + .5)) +xt ;
                if (pym > 0.0 ) {
                      whys[2] = (int) (fact*(-pym -50.)) +yt ;
                }
                else {
                      whys[2] = (int) (fact*(-pym +50.)) +yt ;
                }
                off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                off1Gg.drawLine(exes[1],whys[1],exes[2],whys[2]) ;
                if (pboflag == 3) {    /* smoke trail  MODS  21 JUL 99 */
                  off1Gg.setColor(Color.green) ;
                  exes[1] = (int) (fact*xm[19][1]) + xt ;
                  whys[1] = (int) (fact*(-ym[19][1])) + yt ;
                  for (i=2 ; i<= nptc; ++i) {
                     exes[0] = exes[1] ;
                     whys[0] = whys[1] ;
                     exes[1] = (int) (fact*xm[19][i]) + xt ;
                     whys[1] = (int) (fact*(-ym[19][i])) + yt ;
                     if ((i-antim)/3*3 == (i-antim) ) {
                       off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                     }
                   }
                 }
               }
             }
     // draw the airfoil geometry
             off1Gg.setColor(Color.white) ;
             exes[1] = (int) (fact*(xm[0][npt2])) + xt ;
             whys[1] = (int) (fact*(-ym[0][npt2])) + yt ;
             exes[2] = (int) (fact*(xm[0][npt2])) + xt ;
             whys[2] = (int) (fact*(-ym[0][npt2])) + yt ;
             for (i=1 ; i<= npt2-1; ++i) {
                exes[0] = exes[1] ;
                whys[0] = whys[1] ;
                exes[1] = (int) (fact*(xm[0][npt2-i])) + xt ;
                whys[1] = (int) (fact*(-ym[0][npt2-i])) + yt ;
                exes[3] = exes[2] ;
                whys[3] = whys[2] ;
                exes[2] = (int) (fact*(xm[0][npt2+i])) + xt ;
                whys[2] = (int) (fact*(-ym[0][npt2+i])) + yt ;
                camx[i] = (exes[1] + exes[2]) / 2 ;
                camy[i] = (whys[1] + whys[2]) / 2 ;
                if (foil == 3) {
                    off1Gg.setColor(Color.yellow) ;
                    off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                }
                else {
                    off1Gg.setColor(Color.white) ;
                    off1Gg.fillPolygon(exes,whys,4) ;
                }
             }
             off1Gg.setColor(Color.green) ;
             off1Gg.drawString("View - Side",10,10) ;
             off1Gg.drawString("Flow",10,175) ;
             off1Gg.drawLine(50,172,80,172) ;
             exes[0] = 80 ;  exes[1] = 80; exes[2] = 90 ;
             whys[0] = 177 ;  whys[1] = 167 ; whys[2] = 172  ;
             off1Gg.fillPolygon(exes,whys,3) ;
   
      // put some info on the geometry
             if (displ == 3) {
                inmax = 1 ;
                for (n=1; n <= nptc; ++n) {
                  if(xm[0][n] > xm[0][inmax]) inmax = n ;
                }
                off1Gg.setColor(Color.green) ;
                exes[0] = (int) (fact*(xm[0][inmax])) + xt ;
                whys[0] = (int) (fact*(-ym[0][inmax])) + yt ;
                off1Gg.drawLine(exes[0],whys[0],exes[0]-250,whys[0]) ;
                off1Gg.drawString("Reference",exes[0]-300,whys[0]) ;
                off1Gg.drawString("Angle",exes[0]+20,whys[0]) ;

                off1Gg.setColor(Color.cyan) ;
                exes[1] = (int) (fact*(xm[0][inmax] -
                      4.0*Math.cos(convdr*alfval)))+xt;
                whys[1] = (int) (fact*(-ym[0][inmax] -
                      4.0*Math.sin(convdr*alfval)))+yt;
                off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                off1Gg.drawString("Chord Line",exes[0]+20,whys[0]+20) ;

                off1Gg.setColor(Color.red) ;
                off1Gg.drawLine(exes[1],whys[1],camx[5],camy[5]) ;
                for (i=7 ; i<= npt2-6; i = i+2) {
                   off1Gg.drawLine(camx[i],camy[i],camx[i+1],camy[i+1]) ;
                }
                off1Gg.drawString("Mean Camber Line",exes[0]-70,whys[1]-10) ;
             }
          }

          g.drawImage(offImg1,0,0,this) ;   
       }
    } // end Viewer
  }  // end lokpnl
 
  class Slvpnl extends Panel {
     Foilvel outerparent ;
     Conpnl conpnl ;
     Outpnl outpnl ;

     Slvpnl (Foilvel target) { 
         outerparent = target ;
         setLayout(new GridLayout(2,1,5,5)) ;

         conpnl  = new Conpnl(outerparent) ;
         outpnl  = new Outpnl(outerparent) ;

         add(outpnl) ;
         add(conpnl) ;
     }

     class Conpnl extends Panel {
        Foilvel outerparent ;
        Choice untch ;
        TextField out ;
   
        Conpnl (Foilvel target) { 
          outerparent = target ;
          setLayout(new GridLayout(7,3,5,5)) ;
 
          untch = new Choice() ;
          untch.addItem("Pounds") ;
          untch.addItem("Newtons");
          untch.setBackground(Color.white) ;
          untch.setForeground(Color.red) ;
          untch.select(0) ;

          out = new TextField("12.5",5) ;
          out.setBackground(Color.black) ;
          out.setForeground(Color.green) ;

          add(new Label("Lift", Label.RIGHT)) ;
          add(out) ;
          add(untch) ;

          add(new Label(" ", Label.CENTER)) ;
          add(new Label(" ", Label.CENTER)) ;
          add(new Label("Select Units", Label.CENTER)) ;

          add(new Label(" ", Label.CENTER)) ;
          add(new Label(" ", Label.CENTER)) ;
          add(new Label(" ", Label.CENTER)) ;

          add(new Label("Red dot", Label.RIGHT)) ;
          add(new Label("shows your", Label.CENTER)) ;
          add(new Label("current value", Label.LEFT)) ;

          add(new Label(" ", Label.CENTER)) ;
          add(new Label(" ", Label.CENTER)) ;
          add(new Label(" ", Label.CENTER)) ;

          add(new Label(" ", Label.CENTER)) ;
          add(new Label(" ", Label.CENTER)) ;
          add(new Label(" ", Label.CENTER)) ;

          add(new Label(" ", Label.CENTER)) ;
          add(new Label(" ", Label.CENTER)) ;
          add(new Label(" ", Label.CENTER)) ;
        }

        public boolean action(Event evt, Object arg) {
          if(evt.target instanceof Choice) {
             lunits  = untch.getSelectedIndex() ;
                   // **** the lunits check MUST come first
             setUnits () ;
             this.loadProb () ;

             return true ;
          }
          else return false ;
        } // Handler
 
        public void loadProb() {   // load the input panels
             int i1,i2,i3,i4,i5 ;
             double v1,v2,v3,v4,v5 ;
             float fl1,fl2,fl3,fl4,fl5 ;
                  //  dimensional
             if (lunits == 0) {
                 l.inppnl.lft.l1.setText("Speed-mph") ;
             }
             if (lunits == 1) {
                 l.inppnl.lft.l1.setText("Speed-km/h") ;
             }
             v4 = vfsd ;
             vmn = 0.0;   vmx= vmaxa ;

             fl4 = (float) v4 ;
      
             l.inppnl.lft.f1.setText(String.valueOf(fl4)) ;
      
             i4 = (int) (((v4 - vmn)/(vmx-vmn))*1000.) ;

             l.inppnl.rght.s1.setValue(i4) ;

             computeFlow() ;
             return ;
        }
     } // Conppnl

     class Outpnl extends Canvas
         implements Runnable{
        Foilvel outerparent ;
        Thread run2 ;
        Point locp,ancp;

        Outpnl (Foilvel target) { 
           setBackground(Color.blue) ;
           run2 = null ;
        }

        public void start() {
           if (run2 == null) {
              run2 = new Thread(this) ;
              run2.start() ;
           }
        }

        public void run() {
          int timer ;
 
          timer = 100 ;
          while (true) {
             try { Thread.sleep(timer); }
             catch (InterruptedException e) {}
             slvpnl.outpnl.repaint() ;
          }
        }

        public void loadPlot() {
          double rad,ang,xc,yc,lftref,clref ;
          double del,spd,awng,ppl,tpl,hpl,angl,thkpl,campl,clpl ;
          int index,ic ;

          lines = 1 ;
          clref =  getClplot(camval,thkval,alfval) ;
          if (Math.abs(clref) <= .001) clref = .001 ;    /* protection */
          lftref = clref * q0 * area/lconv/lconv ;
   
          if (dispp == 0) {    // pressure variation
              npt = npt2 ;
              ntr = 3 ;
              nord = nabs = 1 ;
              for (index = 1; index <= npt; ++ index) {
                  pltx[0][index] = xm[0][npt2-index + 1]/4.0 ;
                  plty[0][index] = plp[npt2-index + 1] ;
                  pltx[1][index] = xm[0][npt2+index - 1]/4.0 ;
                  plty[1][index] = plp[npt2+index - 1] ;
                  pltx[2][index] = xm[0][npt2+index - 1]/4.0 ;
                  plty[2][index] = ps0/2116. * pconv ;
              }
              begx = -.5 ;
              endx = .5 ;
              ntikx = 5 ;
              ntiky = 5 ;
              endy=1.02 * ps0/2116. * pconv ;
              begy=.95 * ps0/2116. * pconv ;
              laby = String.valueOf("Press");
              if (lunits == 0) labyu = String.valueOf("psi");
              if (lunits == 1) labyu = String.valueOf("k-Pa");
              labx = String.valueOf(" X ");
              labxu = String.valueOf("% chord");
          }
          if (dispp == 1) {    // velocity variation
              npt = npt2 ;
              ntr = 3 ;
              nord = 2 ;
              nabs = 1 ;
              for (index = 1; index <= npt; ++ index) {
                  pltx[0][index] = xm[0][npt2-index+1]/4.0 ;
                  plty[0][index] = plv[npt2-index+1];
                  pltx[1][index] = xm[0][npt2+index-1]/4.0 ;
                  plty[1][index] = plv[npt2+index-1] ;
                  pltx[2][index] = xm[0][npt2-index+1]/4.0 ;
                  plty[2][index] = vfsd ;
              }
              begx = -.5 ;
              endx = .5 ;
              ntikx = 5 ;
              ntiky = 6 ;
              begy = 0.0 ;
              endy = 500. ;
              laby = String.valueOf("Vel");
              if (lunits == 0) labyu = String.valueOf("mph");
              if (lunits == 1) labyu = String.valueOf("kmh");
              labx = String.valueOf(" X ");
              labxu = String.valueOf("% chord");
          }
          if (dispp == 2) {    // lift versus angle
              npt = 20 ;
              ntr = 1 ;
              nabs = 2;  nord = 3 ;
              begx=-20.0; endx=20.0; ntikx=5;
              labx = String.valueOf("Angle ");
              labxu = String.valueOf("degrees");
              del = 40.0 / npt ;
              for (ic=1; ic <=npt; ++ic) {
                   angl = -20.0 + (ic-1)*del ;
                   clpl = getClplot(camval,thkval,angl) ;
                   pltx[0][ic] = angl ;
                   plty[0][ic] = fconv*lftref * clpl / clref ;
              }
              ntiky = 5 ;
              laby = String.valueOf("Lift");
              if (lunits == 0) labyu = String.valueOf("lbs");
              if (lunits == 1) labyu = String.valueOf("N");
          }
          if (dispp == 3) {    // lift versus thickness
              npt = 20 ;
              ntr = 1 ;
              nabs = 3;  nord = 3 ;
              begx=0.0; endx=25.0; ntikx=6;
              labx = String.valueOf("Thickness ");
              labxu = String.valueOf("% chord");
              del = 1.0 / npt ;
              for (ic=1; ic <=npt; ++ic) {
                   thkpl = .05 + (ic-1)*del ;
                   clpl = getClplot(camval,thkpl,alfval) ;
                   pltx[0][ic] = thkpl*25. ;
                   plty[0][ic] = fconv*lftref * clpl / clref ;
              }
              ntiky = 5 ;
              laby = String.valueOf("Lift");
              if (lunits == 0) labyu = String.valueOf("lbs");
              if (lunits == 1) labyu = String.valueOf("N");
          }
          if (dispp == 4) {    // lift versus camber
              npt = 20 ;
              ntr = 1 ;
              nabs = 4;  nord = 3 ;
              begx=-12.5; endx=12.5; ntikx=5;
              labx = String.valueOf("Camber ");
              labxu = String.valueOf("% chord");
              del = 2.0 / npt ;
              for (ic=1; ic <=npt; ++ic) {
                  campl = -1.0 + (ic-1)*del ;
                  clpl = getClplot(campl,thkval,alfval) ;
                  pltx[0][ic] = campl*12.5 ;
                  plty[0][ic] = fconv*lftref * clpl / clref ;
              }
              ntiky = 5 ;
              laby = String.valueOf("Lift");
              if (lunits == 0) labyu = String.valueOf("lbs");
              if (lunits == 1) labyu = String.valueOf("N");
          }
          if (dispp == 5) {    // lift versus speed
              npt = 20 ;
              ntr = 1 ;
              nabs = 5;  nord = 3 ;
              begx=0.0; 
              labx = String.valueOf("Speed ");
              if (lunits == 0) {
                  labxu = String.valueOf("mph");
                  endx=300.0; ntikx=7;
              }
              if (lunits == 1) {
                  labxu = String.valueOf("kmh");
                  endx=400.0; ntikx=6;
              }
              del = vmaxa / npt ;
              for (ic=1; ic <=npt; ++ic) {
                  spd = (ic-1)*del ;
                  pltx[0][ic] = spd ;
                  plty[0][ic] = fconv*lftref * spd * spd / (vfsd * vfsd) ;
              }
              pltx[1][0] = vfsd ;
              plty[1][0] = lift ;
              begy=0.0;
              ntiky = 5 ;
              laby = String.valueOf("Lift");
              if (lunits == 0) {
                  labyu = String.valueOf("lbs");
                  endy=2000.0; 
              }
              if (lunits == 1) {
                  labyu = String.valueOf("N");
                  endy=8000.0; 
              }
          }
          if (dispp == 6) {    // lift versus altitude
              npt = 20 ;
              ntr = 1 ;
              nabs = 6;  nord = 3 ;
              begx=0.0; endx=50.0; ntikx=6;
              if (lunits == 0) endx = 50.0 ;
              if (lunits == 1) endx = 15.0 ;
              labx = String.valueOf("Altitude");
              if (lunits == 0) labxu = String.valueOf("k-ft");
              if (lunits == 1) labxu = String.valueOf("km");
              del = altmax / npt ;
              for (ic=1; ic <=npt; ++ic) {
                  hpl = (ic-1)*del ;
                  pltx[0][ic] = lconv*hpl/1000. ;
                  tpl = 518.6 ;
                  ppl = 2116. ;
                  if (planet == 0) {
                      if (hpl < 36152.)   {
                            tpl = 518.6 - 3.56 * hpl /1000. ;
                            ppl = 2116. * Math.pow(tpl/518.6, 5.256) ;
                      }
                      else {
                            tpl = 389.98 ;
                            ppl = 2116. * .236 * Math.exp((36000.-hpl)/(53.35*tpl)) ;
                      }
                      plty[0][ic] = fconv*lftref * ppl/(tpl*53.3*32.17) / rho ;
                  }
                  if (planet == 1) {
                      if (hpl <= 22960.) {
                         tpl = 434.02 - .548 * hpl/1000. ;
                         ppl = 14.62 * Math.pow(2.71828,-.00003 * hpl) ;
                      }
                      if (hpl > 22960.) {
                         tpl = 449.36 - 1.217 * hpl/1000. ;
                         ppl = 14.62 * Math.pow(2.71828,-.00003 * hpl) ;
                      }
                      plty[0][ic] = fconv*lftref * ppl/(tpl*1149.) / rho ;
                  }
              }
              ntiky = 5 ;
              laby = String.valueOf("Lift");
              if (lunits == 0) labyu = String.valueOf("lbs");
              if (lunits == 1) labyu = String.valueOf("N");
          }
          if (dispp == 7) {    // lift versus area
              npt = 20 ;
              ntr = 1 ;
              nabs = 7;  nord = 3 ;
              begx=0.0; endx=1000.0; ntikx=6;
              labx = String.valueOf("Area ");
              if (lunits == 0) labxu = String.valueOf("sq ft");
              if (lunits == 1) labxu = String.valueOf("sq m");
              del = armax*lconv*lconv / npt ;
              for (ic=1; ic <=npt; ++ic) {
                   awng = (ic-1)*del ;
                   pltx[0][ic] = awng ;
                   plty[0][ic] = fconv*lftref * awng /area ;
              }
              ntiky = 5 ;
              laby = String.valueOf("Lift");
              if (lunits == 0) labyu = String.valueOf("lbs");
              if (lunits == 1) labyu = String.valueOf("N");
          }
        }

        public double  getClplot (double camb, double thic, double angl) {
           double beta,rball,xc,yc,rc,gamc,lec,tec,lecm,tecm,crdc ;
           double stfact,number ;
   
           xc = 0.0 ;
           yc = camb / 2.0 ;
           rc = thic/4.0 + Math.sqrt( thic*thic/16.0 + yc*yc + 1.0);
           xc = 1.0 - Math.sqrt(rc*rc - yc*yc) ;
           beta = Math.asin(yc/rc)/convdr ;       /* Kutta condition */
           gamc = 2.0*rc*Math.sin((angl+beta)*convdr) ;
           lec = xc - Math.sqrt(rc*rc - yc*yc) ;
           tec = xc + Math.sqrt(rc*rc - yc*yc) ;
           lecm = lec + 1.0/lec ;
           tecm = tec + 1.0/tec ;
           crdc = tecm - lecm ;
                                       // stall model 1
           stfact = 1.0 ;
           if (anflag == 1) {
               if (angl > 10.0 ) {
                  stfact = .5 + .1 * angl - .005 * angl * angl ;
               }
               if (angl < -10.0 ) {
                  stfact = .5 - .1 * angl - .005 * angl * angl ;
               }
           }
   
           number = stfact*gamc*4.0*3.1415926/crdc ;
           return (number) ;
        }

        public void update(Graphics g) {
           slvpnl.outpnl.paint(g) ;
        }
 
        public void paint(Graphics g) {
           int i,j,k,n,index ;
           int xlabel,ylabel,ind,inmax,inmin ;
           int exes[] = new int[8] ;
           int whys[] = new int[8] ;
           double offx,scalex,offy,scaley,waste,incy,incx;
           double xtrans,ytrans,xl,yl;
           double liftab ;
           int camx[] = new int[19] ;
           int camy[] = new int[19] ;
           Color col ;

           off2Gg.setColor(Color.blue) ;
           off2Gg.fillRect(0,0,300,300) ;
           off2Gg.setColor(Color.white) ;
           off2Gg.drawString("Plotter",250,15) ;

           xtrans = 0. ;
           ytrans = 0. ;

           if (ntikx < 2) ntikx = 2 ;     /* protection 13June96 */
           if (ntiky < 2) ntiky = 2 ;
           offx = 0.0 - begx ;
           scalex = 6.0/(endx-begx) ;
           incx = (endx-begx)/(ntikx-1);
           offy = 0.0 - begy ;
           scaley = 4.5/(endy-begy) ;
           incy = (endy-begy)/(ntiky-1) ;

           if (dispp <= 7) {             /*  draw a graph */
                                              /* draw axes */
              off2Gg.setColor(Color.white) ;
              exes[0] = (int) (factp* 0.0 + xtrans) + xtp ;
              whys[0] = (int) (factp* -4.5 + ytrans) + ytp ;
              exes[1] = (int) (factp* 0.0 + xtrans) + xtp ;
              whys[1] = (int) (factp* 0.0 + ytrans) + ytp ;
              exes[2] = (int) (factp* 6.0 + xtrans) + xtp ;
              whys[2] = (int) (factp* 0.0 + ytrans) + ytp ;
              off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
              off2Gg.drawLine(exes[1],whys[1],exes[2],whys[2]) ;

              xlabel = (int) (-90.0 + xtrans) + xtp ;   /*  label y axis */
              ylabel = (int) (factp*-1.5 + ytrans) + ytp ;
              off2Gg.drawString(laby,xlabel,ylabel) ;
              off2Gg.drawString(labyu,xlabel,ylabel+10) ;
                                                    /* add tick values */
              for (ind= 1; ind<= ntiky; ++ind){
                   xlabel = (int) (-50.0+xtrans) + xtp ;
                   yl = begy + (ind-1) * incy ;
                   ylabel = (int) (factp* -scaley*(yl + offy) + ytrans) + ytp ;
                   if (nord >= 2) {
                      off2Gg.drawString(String.valueOf((int) yl),xlabel,ylabel) ;
                   }
                   else {
                      off2Gg.drawString(String.valueOf(filter3(yl)),xlabel,ylabel);
                   }
              }
              xlabel = (int) (factp*3.0 + xtrans) + xtp ;    /* label x axis */
              ylabel = (int) (40.0 + ytrans) + ytp ;
              off2Gg.drawString(labx,xlabel,ylabel-10) ;
              off2Gg.drawString(labxu,xlabel,ylabel) ;
                                                    /* add tick values */
              for (ind= 1; ind<= ntikx; ++ind){
                   ylabel = (int) (15. + ytrans) + ytp ;
                   xl = begx + (ind-1) * incx ;
                   xlabel = (int) (factp*(scalex*(xl + offx) -.05) +xtrans) + xtp ;
                   if (nabs == 1) {
                      off2Gg.drawString(String.valueOf(xl),xlabel,ylabel) ;
                   }
                   if (nabs > 1) {
                      off2Gg.drawString(String.valueOf((int) xl),xlabel,ylabel) ;
                   }
              }
       
              if(lines == 0) {
                 for (i=1; i<=npt; ++i) {
                     xlabel = (int) (factp*scalex*(offx+pltx[0][i])+xtrans) + xtp ;
                     ylabel = (int)(factp*-scaley*(offy+plty[0][i])+ytrans +7.)+ytp;
                     off2Gg.drawString("*",xlabel,ylabel) ;
                 }
              }
              else {
                if (dispp <= 1) {
                 if (anflag == 0 || (anflag == 1 && Math.abs(alfval) < 10.0)) {
                   for (j=0; j<=ntr-1; ++j) {
                      if (j == 0) {
                        off2Gg.setColor(Color.white) ;
                        xlabel = (int) (factp* 6.1 + xtrans) + xtp ;
                        ylabel = (int) (factp* -2.5 + ytrans) + ytp ;
                        off2Gg.drawString("Upper",xlabel,ylabel) ;
                      }
                      if (j == 1) {
                        off2Gg.setColor(Color.yellow) ;
                        xlabel = (int) (factp* 6.1 + xtrans) + xtp ;
                        ylabel = (int) (factp* -1.5 + ytrans) + ytp ;
                        off2Gg.drawString("Lower",xlabel,ylabel) ;
                      }
                      if (j == 2) {
                        off2Gg.setColor(Color.green) ;
                        xlabel = (int) (factp* 2.0 + xtrans) + xtp ;
                        ylabel = (int) (factp* -5.0 + ytrans) + ytp ;
                        off2Gg.drawString("Free Stream",xlabel,ylabel) ;
                      }
                      exes[1] = (int) (factp*scalex*(offx+pltx[j][1])+xtrans) + xtp;
                      whys[1] = (int) (factp*-scaley*(offy+plty[j][1])+ytrans)+ ytp;
                      for (i=1; i<=npt; ++i) {
                        exes[0] = exes[1] ;
                        whys[0] = whys[1] ;
                        exes[1] = (int)(factp*scalex*(offx+pltx[j][i])+xtrans)+xtp;
                        whys[1] = (int)(factp*-scaley*(offy+plty[j][i])+ytrans)+ytp;
                        off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                      }
                   }
                 }
                 if (anflag == 1 && Math.abs(alfval) > 10.0) {
                     off2Gg.setColor(Color.yellow) ;
                     xlabel = (int) (factp* 1.0 + xtrans) + xtp ;
                     ylabel = (int) (factp* -2.0 + ytrans) + ytp ;
                     off2Gg.drawString("Wing is Stalled",xlabel,ylabel) ;
       
                     xlabel = (int) (factp* 1.0 + xtrans) + xtp ;
                     ylabel = (int) (factp* -1.0 + ytrans) + ytp ;
                     off2Gg.drawString("Plot not Available",xlabel,ylabel) ;
                 }
               }
               if (dispp > 1) {
                 off2Gg.setColor(Color.white) ;
                 exes[1] = (int) (factp*scalex*(offx+pltx[0][1])+xtrans) + xtp;
                 whys[1] = (int) (factp*-scaley*(offy+plty[0][1])+ytrans) + ytp;
                 for (i=1; i<=npt; ++i) {
                   exes[0] = exes[1] ;
                   whys[0] = whys[1] ;
                   exes[1] = (int) (factp*scalex*(offx+pltx[0][i])+xtrans) + xtp;
                   whys[1] = (int) (factp*-scaley*(offy+plty[0][i])+ytrans) + ytp;
                   off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                 }
                 xlabel = (int) (factp*scalex*(offx+pltx[1][0])+xtrans) + xtp ;
                 ylabel = (int)(factp*-scaley*(offy+plty[1][0])+ytrans)+ytp -4;
                 off2Gg.setColor(Color.red) ;
//                 off2Gg.drawString("*",xlabel,ylabel) ;
                 off2Gg.fillOval(xlabel,ylabel,5,5) ;
               }
             }
          }
          else {      /*  draw the lift gauge */
              off2Gg.setColor(Color.black) ;
              off2Gg.fillRect(0,100,300,30) ;
            // Thermometer gage
              off2Gg.setColor(Color.white) ;
              if (lftout == 0) off2Gg.drawString("Lift =",70,75) ;
              if (lftout == 1) off2Gg.drawString(" Cl  =",70,75) ;
              off2Gg.setColor(Color.yellow);
              for (index=0 ; index <= 10; index ++) {
                off2Gg.drawLine(7+index*25,100,7+index*25,110) ;
                off2Gg.drawString(String.valueOf(index),5+index*25,125) ;
              }

              if (Math.abs(lift) <= 1.0) {
                 liftab = lift*10.0 ;
                 off2Gg.setColor(Color.cyan);
                 off2Gg.fillRect(0,100,7 + (int) (25*Math.abs(liftab)),10) ;
                 off2Gg.drawString(String.valueOf(filter3(liftab)),110,75) ;
                 off2Gg.drawString(" X 10 ",150,75) ;
                 off2Gg.drawString("-1",180,70) ;
              }
              if (Math.abs(lift) > 1.0 && Math.abs(lift) <= 10.0) {
                 liftab = lift ;
                 off2Gg.setColor(Color.yellow);
                 off2Gg.fillRect(0,100,7 + (int) (25*Math.abs(liftab)),10) ;
                 off2Gg.drawString(String.valueOf(filter3(liftab)),110,75) ;
                 off2Gg.drawString(" X 10 ",150,75) ;
                 off2Gg.drawString("0",180,70) ;
              }
              if (Math.abs(lift) > 10.0 && Math.abs(lift) <=100.0) {
                 liftab = lift/10.0 ;
                 off2Gg.setColor(Color.green);
                 off2Gg.fillRect(0,100,7 + (int) (25*Math.abs(liftab)),10) ;
                 off2Gg.drawString(String.valueOf(filter3(liftab)),110,75) ;
                 off2Gg.drawString(" X 10 ",150,75) ;
                 off2Gg.drawString("1",180,70) ;
              }
              if (Math.abs(lift) > 100.0 && Math.abs(lift) <=1000.0) {
                 liftab = lift/100.0 ;
                 off2Gg.setColor(Color.red);
                 off2Gg.fillRect(0,100,7 + (int) (25*Math.abs(liftab)),10) ;
                 off2Gg.drawString(String.valueOf(filter3(liftab)),110,75) ;
                 off2Gg.drawString(" X 10 ",150,75) ;
                 off2Gg.drawString("2",180,70) ;
              }
              if (Math.abs(lift) > 1000.0 && Math.abs(lift) <=10000.0) {
                 liftab = lift/1000.0 ;
                 off2Gg.setColor(Color.magenta);
                 off2Gg.fillRect(0,100,7 + (int) (25*Math.abs(liftab)),10) ;
                 off2Gg.drawString(String.valueOf(filter3(liftab)),110,75) ;
                 off2Gg.drawString(" X 10 ",150,75) ;
                 off2Gg.drawString("3",180,70) ;
              }
              if (Math.abs(lift) > 10000.0 && Math.abs(lift) <=100000.0) {
                 liftab = lift/10000.0 ;
                 off2Gg.setColor(Color.orange);
                 off2Gg.fillRect(0,100,7 + (int) (25*Math.abs(liftab)),10) ;
                 off2Gg.drawString(String.valueOf(filter3(liftab)),110,75) ;
                 off2Gg.drawString(" X 10 ",150,75) ;
                 off2Gg.drawString("4",180,70) ;
              }
              if (Math.abs(lift) > 100000.0 && Math.abs(lift) <=1000000.0) {
                 liftab = lift/100000.0 ;
                 off2Gg.setColor(Color.white);
                 off2Gg.fillRect(0,100,7 + (int) (25*Math.abs(liftab)),10) ;
                 off2Gg.drawString(String.valueOf(filter3(liftab)),110,75) ;
                 off2Gg.drawString(" X 10 ",150,75) ;
                 off2Gg.drawString("5",180,70) ;
              }
              if (Math.abs(lift) > 1000000.0) {
                 liftab = lift/1000000.0 ;
                 off2Gg.setColor(Color.white);
                 off2Gg.fillRect(0,100,7 + (int) (25*Math.abs(liftab)),10) ;
                 off2Gg.drawString(String.valueOf(filter3(liftab)),110,75) ;
                 off2Gg.drawString(" X 10 ",150,75) ;
                 off2Gg.drawString("6",180,70) ;
              }
          }
    
          g.drawImage(offImg2,0,0,this) ;   
       }
     }  // Outpnl 
  } //Slvpnl 
}
