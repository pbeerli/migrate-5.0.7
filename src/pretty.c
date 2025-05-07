/*
 *  pretty.c
 *  migrate-n
 *
 *  Created by Peter Beerli on 7/25/05.
 *  changed to accomodate libharu (libhpdf) 2.4 August 2023
 *  Copyright 2005-2023 Peter Beerli. All rights reserved.
 
 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject
 to the following conditions:
 
 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
 *
 */
//#ifdef PRETTY
#pragma clang diagnostic ignored "-Wformat-nonliteral"
#include "pretty.h"
#include "data.h"
#include "options.h"
#include "migevents.h"
#include "mutationmodel.h"
#include "menu.h"
#include "priors.h"
#include "bayes.h"
#include "speciate.h"
#include "world.h"
#include <hpdf.h>
#ifdef BEAGLE
#include "beagle.h"
#endif
#include <stdarg.h>

#define LINESTRETCH 15
//#define NOT_PERCENTILES FALSE
//#define PERCENTILES TRUE
#define SQUARE 0
#define DIAMOND 1

#ifndef A4PAPER
#define LETTERADJUST 50
#else
#define LETTERADJUST 0
#endif

//#ifndef HAS_INDIX
#define INDIX(a,b,c) ((a)*(b)+(c))
//#define HAS_INDIX
//#endif
extern int myID;	  /*!< myID=0 for single-cpu program and master in MPI program, worker nodes myID > 0 */
extern time_t startseconds;
HPDF_Doc doc;
HPDF_Page page;
//HPDF_Page canvas; // new it seems that page and canvas are the same
HPDF_Font helv;
HPDF_Font helvob;
HPDF_Font symb;
HPDF_Font cour;
int page_counter;
char pdf_pagetitle[LINESIZE+1];
char pdf_time[LINESIZE+1];
double page_height;
double left_margin;
extern int numcpu;
static HPDF_Page firstpage;
double timestampy;
void pdf_printf_right(double x, double y, char string[],...);
int pdf_print_header(char *title);
void pdf_print_contents_at(double x, double y, char *title);
void pdf_draw_line(double xs, double ys, double xe, double ye);
void pdf_migrate_logo(double x, double y, double stretch);
void pdf_migrate_logo_lines(double x, double y, double stretch);
void pdf_print_time(double startx, double *pageheight, char text[]);
boolean  pdf_advance(double *page_height);
void pdf_print_seqdata (double margin, world_fmt * world, data_fmt * data, option_fmt * options);
void pdf_printf_next(double x, double *y, char string[], ...);
void pdf_printf_ralign(double rx, double y, char string[], ...);
void pdf_printf(double x, double y, char align, char string[], ...);
void pdf_print_tableline(double *width, char *fmt, ...);
void pdf_print_result_param (double *lx,  MYREAL *param, long numpop, long pop,
                             boolean usem);
void pdf_table(int cols, int rows, char **buffer, char *position, int col_overflow, double separator);

long nice_element(MYREAL param, char *element, MYREAL lower, MYREAL mid, MYREAL upper,
                  int low_mid_digits, int mid_upper_digits, char delimiter);
void pdf_linedotplot(long n, double *x, double *y, double dotcolor[], double linecolor[], MYREAL linethickness, double width, double height, boolean has_dots, boolean has_line);

void pdf_print_symbol(double lx, int fontsize, char pos, char *symbolstring);

void pdf_print_citation(char type[],world_fmt *world);
void pdf_print_section_title(double *page_width, double * page_height, char *title);
double pdf_print_line_element(double lx, double ly, double offset, char *title);
double pdf_print_line_element2(double lx, double ly, double offset, double value, int fmt1, int fmt2);
void symbol_Theta(double lx, double ly, int size, long subscript);
void symbol_Growth(double lx, double ly, int size, long subscript);
void symbol_mlf(double lx, double ly, int size, long subscript);
double quantiler(double *values, double prob1, double prob2, long range1, long range2);
void symbol_R(double lx, double ly, int size, long subscript);
void symbol_M(double lx, double ly, int size, long subscript1, long subscript2, boolean usem);
void symbol_D(double lx, double ly, int size, long subscript1, long subscript2);
void symbol_S(double lx, double ly, int size, long subscript1, long subscript2);
void symbol_Hexp(double lx, double ly, int size);
void draw_rect(HPDF_Page page, double x, double y, const char* label);
double pdf_page_advance_or_not(double *mypage_height, double need_pixels);
boolean pdf_advance_half(double *mypage_height);
void   pdf_draw_tick(double xs, double ys, int orientation, double ticklength, double value, int digits);
double prettytick(double mi, double ma, long ticks,
                 double *nmi, double *nma, double * newdelta, int * digits);
void  pdf_create_axes(int up, double xmin, double xmax, long ticks, double lx, double ly, double length);
void pdf_print_dot(double xs, double ys, double width, int shape, double color[]);
void findmoments(double *vals, long n, double *p01, long *pos01, double *p99, long *pos99, double *mm, long *posmm);
void findoutliers(double *vals, long n,
                  double *p01, long *pos01,
                  double *p99, long *pos99,
                  double *mm, long *posmm);
void findminmax_freq_only(double *vals, long n, double *p00, long *pos00, double *p100, long *pos100);
void pdf_print_line_species(char *title, double lx, double ly, double offset, long frompop, long topop);
void  pdf_print_line_mig(char *migtitle, double lx, double page_height, double offset, long frompop, long topop);
void pdf_print_bayestableheader(double mypage_height, double left_margin, double right_margin, double *offset);
void pdf_print_line_theta(double lx, double ly, double offset, long j);
void pdf_print_line_growth(double lx, double ly, double offset, long j);
void pdf_print_line_mlf(double lx, double ly, double offset, long j);
void pdf_print_line_rate(double lx, double ly, double offset, long j, long exponent);
void pdf_fill_stroke(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
void pdf_title(char *title, double page_width);
void find_posterior_min_max(double *minval, double *maxval, long start, long stop, bayes_fmt *bayes, long locus);
void pdf_pretty_histogram(long pa, long rpa,long numbins, double * results, long stride, char * set50,
                          char * set95, long * bins, double delta, MYREAL * mini, MYREAL * maxi,
                          double lx, double ly, world_fmt *world, double themin, double themax);
void pdf_advance_hist(long numpop, long *z, double *lx, double *ly);
void pdf_putc(double *x, double *y, double leftborder, double rightborder, char message);
void pdf_printf_cell(double *x, double *y, double width, char string[], ...);
void    pdf_print_connection_table (world_fmt * world, option_fmt *options, data_fmt * data);
void    pdf_print_connection_table (world_fmt * world, option_fmt *options, data_fmt * data);
void pdf_printf_right_next(double x, double *y, char string[],...);
void    pdf_print_distance_table (world_fmt * world, option_fmt * options, data_fmt * data);
void pdf_print_ratetbl (world_fmt * world, option_fmt * options, long locus, char header);
void pdf_print_param_order(world_fmt *world);
int count_elements(char *fmt);
void pdf_print_allelelegend(double *column_width, long loci);
void pdf_print_dataheader (boolean first, char *title, double *column_width, long pop, world_fmt * world,
			   option_fmt * options, data_fmt * data);
void pdf_print_alleledata (double margin, world_fmt * world, data_fmt * data, option_fmt * options);
void pdf_print_sequence(double right_margin, data_fmt *data, long locus, long pop, long ind);
void format_helper(MYREAL * matrix, long len, int *fmt1, int *fmt2);
void pdf_print_matrix_line(long whichline, double width, long cols,long rows,MYREAL * matrix,boolean ismig);
void pdf_print_result_header (double *lx, char *titletext, world_fmt * world);
void pdf_print_popstring(double *lx, long pop, world_fmt *world, option_fmt *options, data_fmt *data);
void pdf_print_replicate(double lx, world_fmt *world, long maxrep, long rep, long locus);
void pdf_print_result_population (double *lx, long pop, world_fmt * world,
				  option_fmt * options, data_fmt * data);
void pdf_print_comment(double lx, double *ly, char *this_text);
void    pdf_print_correlation_table (world_fmt * world, option_fmt * options, data_fmt * data);
void method_set(double lx, char method);
void pdf_table_footnote(double lx, long failed, boolean percentiles);
void  translate_buffer_table(long cols, long rows, char **thebuffer, char **header, char ***elements);
void  extract_column_buffer_table(long col, long cols, long rows, char **thebuffer, double *x);
void  find_col_width(int cols, int rows, char ***elements, char **header, double *col_widths);
double align_column(char position, double cw, double col_leftmargin);
void  define_col_start(double cols, double * col_widths, int col_overflow, char *position, double page_width, double separator, double *col_starts);
void  pdf_print_table_header(char * position, int cols, double *col_starts,char **header);
void pdf_print_symbol_no(char value, char* symbolstring, double lx, char pos);
long pdf_table2(int cols, int rows, char **header, char **header2, char ***elements, char *position, int col_overflow, double separator);
void pdf_print_endline(void);
void findminmax(double *vals, const long n, double *min, double *max);
//##errorhandler

#include <setjmp.h>

static jmp_buf env;

#ifdef HPDF_DLL
void  __stdcall
#else
void
#endif
error_handler (HPDF_STATUS   error_no,
               HPDF_STATUS   detail_no,
               void         *user_data)  __attribute__((noreturn));

#ifdef HPDF_DLL
void  __stdcall
#else
void
#endif
error_handler (HPDF_STATUS   error_no,
               HPDF_STATUS   detail_no,
               void         *user_data)
{
  (void) user_data;
  printf ("ERROR: error_no=%04X, detail_no=%u\n", (HPDF_UINT)error_no,
	  (HPDF_UINT)detail_no);
  longjmp(env, 1);
}

///
/// returns the quantile at prob1 or prob2, if the value at prob2 is smaller than
/// the value at prob1 then the value at prob1 is return otherwise the value at prob2.
/// We assume range1< range2
double quantiler(double *values, double prob1, double prob2, long range1, long range2)
{
    double *temp;
    MYREAL a,b;
    if (range2==0)
        return 0.0;
    temp = (double *) mycalloc((range2+1),sizeof(double));
    memcpy(temp, values,sizeof(double)*(size_t) range2);
    qsort ((void *) temp, (size_t) range2, sizeof (double), doublecmp);
    a = temp[((long) (range1 * prob1))];
    b = temp[((long) ((range2-1) * prob2))];
    myfree(temp);
    return (double) (a<b ? b : a);
}

void symbol_Theta(double lx, double ly, int size, long subscript)
{
    char *thetatitle = "Q";
    char tempstring[100];
    int subsize =  (int) (0.75 * size);
    double lxdelta = lx + (double) 0.8 * size;
    double lydelta = ly - (double) 0.5 * subsize;
    if(subscript >= 0)
        sprintf(tempstring,"%li",subscript);
    else
        sprintf(tempstring," ");
    HPDF_Page_SetFontAndSize (page, symb, (float) size);
    pdf_print_contents_at(lx,ly,thetatitle);
    HPDF_Page_SetFontAndSize (page, symb, (float) subsize);
    pdf_print_contents_at(lxdelta,lydelta,tempstring);
    // backt to default
    HPDF_Page_SetFontAndSize (page, helv, 10);
}

void symbol_Growth(double lx, double ly, int size, long subscript)
{
    char *growthtitle = "g";
    char tempstring[100];
    int subsize =  (int) (0.75 * size);
    double lxdelta = lx + (double) 0.8 * size;
    double lydelta = ly - (double) 0.5 * subsize;
    if(subscript >= 0)
        sprintf(tempstring,"%li",subscript);
    else
        sprintf(tempstring," ");
    HPDF_Page_SetFontAndSize (page, helv, (float) size);
    pdf_print_contents_at(lx,ly,growthtitle);
    HPDF_Page_SetFontAndSize (page, symb, (float) subsize);
    pdf_print_contents_at(lxdelta,lydelta,tempstring);
    // backt to default
    HPDF_Page_SetFontAndSize (page, helv, 10);
}
void symbol_mlf(double lx, double ly, int size, long subscript)
{
    char *growthtitle = "a";
    char tempstring[100];
    int subsize =  (int) (0.75 * size);
    double lxdelta = lx + (double) 0.8 * size;
    double lydelta = ly - (double) 0.5 * subsize;
    if(subscript >= 0)
        sprintf(tempstring,"%li",subscript);
    else
        sprintf(tempstring," ");
    HPDF_Page_SetFontAndSize (page, symb, (float) size);
    pdf_print_contents_at(lx,ly,growthtitle);
    HPDF_Page_SetFontAndSize (page, symb, (float) subsize);
    pdf_print_contents_at(lxdelta,lydelta,tempstring);
    // backt to default
    HPDF_Page_SetFontAndSize (page, helv, 10);
}

void symbol_R(double lx, double ly, int size, long subscript)
{
    char *thetatitle = "m";
    char tempstring[100];
    //int subsize =  (int) (0.75 * size);
    //double lxdelta = lx + (double) 0.8 * size;
    //double lydelta = ly - (double) 0.5 * subsize;
    if(subscript > -1)
        sprintf(tempstring,"%li",subscript);
    else
    {
        if(subscript < -1)
            sprintf(tempstring,"combined");
        else
            sprintf(tempstring," ");
    }
    HPDF_Page_SetFontAndSize (page, symb, (float) size);
    pdf_print_contents_at(lx,ly,thetatitle);

    HPDF_Page_SetFontAndSize (page, symb, 10.0);
    sprintf(tempstring,"[10");
    pdf_print_contents_at(lx+13,ly,tempstring); //print the scale should look like this [x10-5]
    HPDF_Page_SetFontAndSize (page, symb, 8.0);
    sprintf(tempstring,"%li",subscript);
    pdf_print_contents_at(lx+26,ly+4,tempstring); // print the exponent as superscript
    HPDF_Page_SetFontAndSize (page, symb, 10.0);
    sprintf(tempstring,"]");
    pdf_print_contents_at(lx+35,ly,tempstring);


    //    if(subscript < -1)
    //    HPDF_Page_SetFontAndSize (page, helv, subsize);
    //else
    //    HPDF_Page_SetFontAndSize (page, symb, subsize);
    //pdf_print_contents_at(lxdelta,lydelta,tempstring);
    // backt to default
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
}

void symbol_M(double lx, double ly, int size, long subscript1, long subscript2, boolean usem)
{
    char *migt1 = "M";
    char *migt2 = "xNm";
    char *migtitle = usem ? migt1 : migt2;
    int msub = usem ?  size/2 : size+ size/2;
    char tempstring[100];
    int subsize =  (int) (0.75 * size);
    double lxdelta = lx + (double) 0.8 * size + msub;
    double lydelta = ly - (double) 0.5 * subsize;
    sprintf(tempstring,"%li->%li",subscript1,subscript2);
    HPDF_Page_SetFontAndSize (page, helv, (float) size);
    pdf_print_contents_at(lx,ly,migtitle);
    HPDF_Page_SetFontAndSize (page, symb, (float) subsize);
    pdf_print_contents_at(lxdelta,lydelta,tempstring);
    // backt to default
    HPDF_Page_SetFontAndSize (page, helv, 10);
}

void symbol_D(double lx, double ly, int size, long subscript1, long subscript2)
{
    char *migtitle ="D";
    int msub = size/2;
    char tempstring[100];
    int subsize =  (int) (0.75 * size);
    double lxdelta = lx + (double) 0.8 * size + msub;
    double lydelta = ly - (double) 0.5 * subsize;
    sprintf(tempstring,"%li->%li",subscript1,subscript2);
    HPDF_Page_SetFontAndSize (page, symb, size);
    pdf_print_contents_at(lx,ly,migtitle);
    HPDF_Page_SetFontAndSize (page, symb, subsize);
    pdf_print_contents_at(lxdelta,lydelta,tempstring);
    // backt to default
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
}

void symbol_S(double lx, double ly, int size, long subscript1, long subscript2)
{
    char *migtitle ="s";
    int msub = size/2;
    char tempstring[100];
    int subsize =  (int) (0.75 * size);
    double lxdelta = lx + (double) 0.8 * size + msub;
    double lydelta = ly - (double) 0.5 * subsize;
    sprintf(tempstring,"%li->%li",subscript1,subscript2);
    HPDF_Page_SetFontAndSize (page, symb, (float) size);
    pdf_print_contents_at(lx,ly,migtitle);
    HPDF_Page_SetFontAndSize (page, symb, (float) subsize);
    pdf_print_contents_at(lxdelta,lydelta,tempstring);
    // backt to default
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
}

void symbol_Hexp(double lx, double ly, int size)
{
    char *title = "H";
    char tempstring[100];
    int subsize =  (int) (0.75 * size);
    double lxdelta = lx + (double) 0.8 * size;
    double lydelta = ly - (double) 0.5 * subsize;
    sprintf(tempstring,"exp");
    HPDF_Page_SetFontAndSize (page, helv, (float) size);
    pdf_print_contents_at(lx,ly,title);
    HPDF_Page_SetFontAndSize (page, helv, (float) subsize);
    pdf_print_contents_at(lxdelta,lydelta,tempstring);
    // backt to default
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
}

///
/// Draw a simple line from x0,y0 to x0/y0, this function does not change
/// thickness or line type
void pdf_draw_line(double xs, double ys, double xe, double ye)
{
  HPDF_Page_MoveTo(page, (float) xs, (float) ys);
  HPDF_Page_LineTo(page, (float) xe, (float) ye);
    HPDF_Page_Stroke(page);
}

///
/// draw a rectangle [from haru examples]
void draw_rect(HPDF_Page mypage, double x, double y, const char* label)
{
    HPDF_Page_SetFontAndSize (mypage, helv, 10.0);
    HPDF_Page_BeginText(mypage);
    HPDF_Page_TextOut(mypage, (float) x, (float) (y - 10.), label);
    HPDF_Page_EndText(mypage);
    
    HPDF_Page_Rectangle(mypage, (float) x, (float) (y - 40.), 220., 25.);
}

///
/// when there are fewer than need_pixels create a new page and
/// and push page_height down so that the need_pixel object will fit.
double pdf_page_advance_or_not(double *mypage_height, double need_pixels)
{
    boolean new_page = FALSE;
    if(*mypage_height < need_pixels)
    {
        pdf_new_page("");
        *mypage_height = (double) (HPDF_Page_GetHeight(page));
        *mypage_height -= 55. + LINESTRETCH;
        new_page = TRUE;
    }
    return new_page;
}

///
/// advance a line and check whether we are at the end of the page, if yes then add a new page
boolean pdf_advance(double *mypage_height)
{
    *mypage_height -= LINESTRETCH;
    return (boolean) (pdf_page_advance_or_not(mypage_height, 55.));
}
///
/// advance a half-line and check whether we are at the end of the page, if yes then add a new page
boolean pdf_advance_half(double *mypage_height)
{
    *mypage_height -= LINESTRETCH/2.;
    return (boolean) (pdf_page_advance_or_not(mypage_height, 55.));
}

#define HORIZONTAL 0
#define VERTICAL   1
///
/// plot a single tick with label at location x/y.
// \param xs  x coordinate in true paper coordinates (points)
// \param xy  y coordinate in true paper coordinates (points)
/// \param orientation either HORIZONTAL or VERTICAL, on HORIZONTAL=0 is checked
/// \param ticklength length of a the tick to draw
/// \param value this is a double value to put into the label
/// \param digits number of decimal digits
void   pdf_draw_tick(double xs, double ys, int orientation, double ticklength, double value, int digits)
{
    double w;
    double h;
    double tl = ticklength;
    char *title;
    title = (char *) mycalloc(100,sizeof(char));
    sprintf(title,"%.*f",digits, value);
    w = (double) HPDF_Page_TextWidth(page, title);
    h = (double) HPDF_Page_GetCurrentFontSize(page);
    if(orientation==HORIZONTAL)
    {
        pdf_print_contents_at(xs-w/2, ys-tl-h, title);
        pdf_draw_line(xs,ys, xs , ys - tl);
    }
    else
    {
        pdf_print_contents_at(xs-tl-w-2, ys-h/2+2, title);
        pdf_draw_line(xs,ys, xs-tl , ys);
    }
    myfree(title);
}

///
/// find pretty values for tick labels
/// using heuristics
double prettytick(double mi, double ma, long ticks,
                 double *nmi, double *nma, double * newdelta, int * digits)
{
  (void) nma;
    double newvalue = 0;
    double dist = (ma - mi)/ticks;
    double multiplier[5];
    long  lmi = (long) (floor(log10(dist)));
    double scale  = (double) pow(10.,lmi);
    long i=0;
    multiplier[0] = 1.;
    multiplier[1] = 2.;
    multiplier[2] = 2.5;
    multiplier[3] = 5.;
    multiplier[4] = 10.;
    if(lmi<3)
      *digits = (int) labs(lmi);
    else
        *digits = 0;
    dist /= scale;
    if(dist < 1.)
        dist *= 10.;
    else
        dist /= 10.;
    if(dist >= 10.)
        dist /= 10.;
    else
        dist *= 10.;
    i = 0;
    while(i < 5 && dist > multiplier[i])
    {
        i++;
    }
    *newdelta = scale * multiplier[i];
    *nmi = (double) ceil(mi/(*newdelta) - 0.05) * (*newdelta);
    
    return newvalue;
}

///
/// create a vertical or horizontal axis
/// \param up is either HORIZONTAL=1 or VERTICAL
/// \param xmin  real start value
/// \param xmax  real end value
/// \param ticks number of ticks
/// \param lx    x paper coordinate to start axis
/// \param ly    y paper coordinate to start axis
/// \param length  length on paper in paper coordinate system
void  pdf_create_axes(int up, double xmin, double xmax, long ticks, double lx, double ly, double length)
{
    long i;
    int digits = 3;
    double value = 0.0;
    double plotdelta = length/ticks;
    double realdelta;
    double newxmin;
    double newxmax;
    double newrealdelta;
    double newlx;
    double newly;
    // pretty printing of labels
    // find pretty tickmarks
    if (fabs(xmin-xmax)<DBL_EPSILON)
    {
        if (xmin!=0.0)
        {
            xmax=xmin + 0.1*xmin;
        }
        else
        {
            xmax = 0.1;
        }
    }
    realdelta = (xmax - xmin)/ticks;
    prettytick(xmin,xmax, ticks, &newxmin, &newxmax, &newrealdelta,&digits);
    //
    if(up!=HORIZONTAL)
    {
        //xcode newlx = lx;
        newly = ly + (newxmin-xmin)/(xmax-xmin) * length;
        plotdelta *= newrealdelta/realdelta;
        pdf_draw_line(lx,ly, lx ,ly + length);
        for(i=0; i<ticks; i++)
        {
            value = newxmin+i*newrealdelta;
            if(value>xmax)
                break;
            pdf_draw_tick(lx, newly+i*plotdelta, VERTICAL, 3, value, digits);
        }
    }
    else
    {
        pdf_draw_line(lx,ly, lx + length, ly);
        newly = ly;
        newlx = lx + (newxmin-xmin)/(xmax-xmin) * length;;
        plotdelta *= newrealdelta/realdelta;
        
        for(i=0; i<ticks; i++)
        {
            value = newxmin+i*newrealdelta;
            if(value>xmax)
                break;
            pdf_draw_tick(newlx + i*plotdelta,
                          newly, HORIZONTAL, 3, value, digits);
        }
    }
}

///
/// plot a dot of the shape i (currently i = square or diamond) at position
/// xs and ys in RGB color is a double vector of 3 values
void pdf_print_dot(double xs, double ys, double width, int shape, double color[])
{
  float x = (float) xs;
  float y = (float) ys; 
  float w = (float) (width / 2.);
  float red = (float) color[0];
  float green = (float) color[1];
  float blue = (float) color[2];
  HPDF_Page_SetRGBFill(page, red, green, blue);
  switch(shape)
    {
    case SQUARE:
      HPDF_Page_MoveTo(page, x-w, y-w);
      HPDF_Page_LineTo(page, x+w, y-w);
      HPDF_Page_LineTo(page, x+w, y+w);
      HPDF_Page_LineTo(page, x-w, y+w);
      HPDF_Page_LineTo(page, x-w, y-w);
      break;
    case DIAMOND:
    default:
            HPDF_Page_MoveTo(page, x, y + w);
            HPDF_Page_LineTo(page, x+w, y);
            HPDF_Page_LineTo(page, x, y-w);
            HPDF_Page_LineTo(page, x-w, y);
            HPDF_Page_LineTo(page, x, y + w);
            break;
    }
    HPDF_Page_Fill(page);
    HPDF_Page_SetRGBFill(page, 0, 0, 0);
}

void findmoments(double *vals, long n, double *p01, long *pos01, double *p99, long *pos99, double *mm, long *posmm)
{
    long i;
    double val;
    double total = 0.;
    *mm = - (double) (LONG_MAX);
    *p01 = 0.;
    *p99 = 0.;
    *pos01 = 0;
    *pos99 = 0;
    *posmm = 0;
    for(i=0;i<n; i++)
    {
        val = (double) vals[i];
        total += val;
        if( *mm < val)
        {
            *mm = val;
            *posmm = i;
        }
    }
    val = 0.;
    for(i=0;i<n; i++)
    {
        val += (double) vals[i]/total;
        //printf("findmoments:val:: %li %f %f %f\n",i,val,vals[i],total);
        if(val < 0.01)
        {
            //	    if(*p01  < (double) vals[i])
            //  {
            *p01 = (double) vals[i];
            *pos01 = i;
            //  }
        }
        if(val <= 0.99)
        {
            //	    if(*p99 < (double) vals[i])
            //  {
            *p99 = (double) vals[i];
            *pos99 = i;
            //   }
        }
    }
}

///
/// finds outliers in vals and retunrs some specified quantiles and the position of these
/// values in the vals vector.
void findoutliers(double *vals, long n,
                  double *p01, long *pos01,
                  double *p99, long *pos99,
                  double *mm, long *posmm)
{
    long i;
    double * temp;
    if(n==0)
    {
        warning("findoutliers(): The number of values to check for the histogram is zero\n");
        return;
    }
    temp = (double *) mycalloc((n+1),sizeof(double));
    memcpy(temp,vals,(size_t)n*sizeof(double));
    qsort(temp, (size_t) n, sizeof(double), doublecmp);
    
    *mm = (double) temp[n-1];
    *p01 = (double) temp[(long)(n * 0.01)];
    *p99 = (double)temp[(long)(n * 0.99)];
    *pos01 = 0;
    *pos99 = 0;
    *posmm = 0;
    for(i=0;i<n; i++)
    {
        if(vals[i] <= *p01)
        {
            *pos01 = i;
        }
        if(vals[i] < *p99)
        {
            *pos99 = i;
        }
        
        if(vals[i] < *mm)
        {
            *posmm = i;
        }
        
    }
    myfree(temp);
}

///
/// only for fractions 0..1
void findminmax_freq_only(double *vals, long n, double *p00, long *pos00, double *p100, long *pos100)
{
    long i;
    double val;
    double total = 0.;
    *p00 = 0.;
    *p100 = 0.;
    *pos00 = 0;
    *pos100 = 0;
    for(i=0;i<n; i++)
    {
        val = (double) vals[i];
        total += val;
    }
    val = 0.;
    for(i=0;i<n; i++)
    {
        val += (double) vals[i]/total;
        if(val < 0.000001)
        {
            *p00 = (double) vals[i];
            *pos00 = i;
        }
        if(val < 0.999999)
        {
            *p100 = (double) vals[i];
            *pos100 = i;
        }
    }
}

///
/// create a histogram from bincounts at a specific location and width and height
/// If binmax is -9999 then the maximum is reset to the 99% percentile.
/// If binmax is -999 then the maximum is reset to the 100% percentile
/// if nofreq is TRUE then the y axes is treated literally and not scaled
void pdf_histogram(double *binvals, char *set50, char *set95, long bins, double bindelta, double binmin, double binmax, double lx, double ly, double width, double height, boolean nofreq, MYREAL *priors)
{
    long i;
    //double total=0.0;
    double binvalsmax=0.0;
    //double priorsmax=0.0f;
    long binvalspos;
    double p99=0.0;
    long pos99=0;
    double p100=0.0;
    long pos100=0;
    double p00=0.0;
    long pos00=0;
    double p01=0.0;
    long pos01=0;
    double x = 0.;
    double delta = width / bins;
    long numbins = bins;
    
    double red[3]={0.99,0.,0.};
    //double sumval=0.0;
    //double sumprior=0.0;
#ifdef DEBUG
    //    double ss1=0.0;
    //    double ss2=0.0;
#endif 
    if(nofreq)
    {
        findoutliers(binvals,bins,&p01,&pos01, &p99,
                     &pos99, &binvalsmax, &binvalspos);
    }
    else
        findmoments(binvals,bins,&p01,&pos01, &p99,
                    &pos99, &binvalsmax, &binvalspos);
    //fprintf(stdout,"@#@#@ p01=%f, pos01=%li, p99=%f, pos99=%li, binvalsmax=%f, binvalspos=%li\n",p01,pos01, p99,
    //	    pos99, binvalsmax, binvalspos);
    //total = 0.;
    //if(priors == NULL)
    //{
    //    for(i=0;i<bins;i++)
    //    {
    //        total += (double) binvals[i];
    //    }
    //}
//    else
//    {
//        priorsmax=-HUGE;
//        for(i=0;i<bins;i++)
//        {
//            if (priors[i]>priorsmax)
//                priorsmax = priors[i];
//            total += (double) binvals[i];
//        }
//        priorsmax = priorsmax/binvalsmax;
//    }
    if(binmax < -9000)
    {
        binmax = binmin + pos99 * bindelta;
        delta = width/pos99;
        numbins = pos99;
    }
    else
    {
        if(binmax < -900)
        {
            findminmax_freq_only(binvals,bins,&p00,&pos00, &p100,
                                 &pos100);
            binmax = binmin + pos100 * bindelta;
            delta = width/pos100;
            numbins = pos100;
        }
    }
    if(nofreq)
    {
        if (p99==0.0)
            p99=binmax;
        pdf_create_axes(VERTICAL, 0.0, p99, 5, lx, ly, height);
    }
    else
        pdf_create_axes(VERTICAL, 0, binvalsmax/*(total)**@#@#@#@***/, 5, lx, ly, height);
    pdf_create_axes(HORIZONTAL, binmin , binmax, 5, lx, ly, width);
    HPDF_Page_SetLineWidth(page, (float) delta);
    
    x = (double) (delta/2.0);
    for(i=0;i<numbins;i++)
    {
        if(set50[i] == '1')
            HPDF_Page_SetRGBStroke(page, 0.1f, 0.1f, 0.1f);
        else
        {
            if(set95[i] == '1')
                HPDF_Page_SetRGBStroke(page, 0.4f, 0.4f, 0.4f);
            else
                HPDF_Page_SetRGBStroke(page, 0.8f, 0.8f, 0.8f);
        }
        if(nofreq)
        {
            if(binvals[i]/p99 > 1.0)
            {
                pdf_draw_line(lx+x,ly+0,lx+x,ly + height);
                pdf_print_dot(lx+x, ly + height - 3.0, 3, SQUARE, red);
            }
            else
                pdf_draw_line(lx+x,ly+0,lx+x,ly + binvals[i]/p99 * height);
        }
        else
            pdf_draw_line(lx+x,ly+0,lx+x,ly + binvals[i]/binvalsmax * height);
#ifdef DEBUG
	//fprintf(stderr,"stderr> %li: %f histheight=%f val=%f / %f\n",i, lx+x,ly+ binvals[i]/binvalsmax * height, binvals[i],binvalsmax);       
	//fprintf(stderr,"stderr> %li: %f histheight=%f val=%f / %f\n",i, lx+x,ly+ binvals[i]/binvalsmax * height, binvals[i],binvalsmax);
#endif
        if(priors !=NULL)
        {
            pdf_print_dot(lx+x, ly + (((double) priors[i]/binvalsmax) * height) , 1, SQUARE, red);
        }
        x += delta;
    }
    //printf ("###########>>>>>>>> %f %f\n", sumval, sumprior);
    HPDF_Page_SetRGBStroke(page, 0.f, 0.f, 0.f);
    HPDF_Page_SetLineWidth(page, 1);
}

///
/// create a histogram from bincounts at a specific location and width and height using *std as +-std dev
/// If binmax is -9999 then the maximum y-value printed is reset to the 99% percentile.
/// If binmax is -999 then the maximum y-value printed is reset to the 100% percentile
/// if nofreq is TRUE then the y axes is treated literally and are not scaled
/// world is needed for the skyline plots so that we can plot the dots for
/// the data points onto the graph.
void pdf_histogram_plus(double *binvals, MYREAL *std,
                        char *set50, char *set95,
                        long bins, double bindelta, double binmin, double binmax,
                        double lx, double ly, double width, double height,
                        MYREAL valmax, boolean nofreq, world_fmt * world, double * confidence, long topop)
{
    long locus;
    long ind;
    long pop;
    double ratio;
    double xcoord;
    double averagerate=0.;
    long i;
    //long nbmin;
    double total=0.0;
    double binvalsmax=0.0;
    long binvalspos;
    double p99=0.0;
    long pos99= -1;
    double p100=0.0;
    long pos100= -1;
    double p00=0.0;
    long pos00= -1;
    double p01=0.0;
    long pos01= -1;
    double x = 0.0;
    double delta = width / bins;
    long numbins = bins;
    //long grouping = 0;
    float r;
    double vold=0.;
    double xold;
    double *upper;
    double *lower;
    MYREAL l;
    MYREAL v;
    MYREAL u;
    //    long k;
    double blue[3]={0.0,0.0,0.99};
    double red[3]={0.99,0.0,0.0};
    
    upper = (double*) mycalloc(numbins,sizeof(double));
    lower = (double*) mycalloc(numbins,sizeof(double));
    // setting lower and upper to the mean values;
    memcpy(upper, binvals, (size_t) numbins * sizeof(double));
    memcpy(lower, binvals, (size_t) numbins * sizeof(double));
    
    for(i=0; i< numbins; i++)
    {
        upper[i] += fabs(std[i]);
        lower[i] -= fabs(std[i]);
        if(lower[i] < 0.0)
            lower[i] = 0.0;
    }
    if(nofreq)
    {
        findoutliers(binvals,bins,&p01,&pos01, &p99,
                     &pos99, &binvalsmax, &binvalspos);
    }
    else
        findmoments(upper,bins,&p01,&pos01, &p99,
                    &pos99, &binvalsmax, &binvalspos);
    // fprintf(stdout,"p01=%f, pos01=%li, p99=%f, pos99=%li, binvalsmax=%f, binvalspos=%li\n",p01,pos01, p99,
    //	    pos99, binvalsmax, binvalspos);
    total = 0.;
    for(i=0;i<bins;i++)
    {
        total += (double) upper[i];
    }
    
    if(binmax < -9000)
    {
        //find the 99% percentile of the upper limit to show the whole skyline
        //the closeup should not use this setting because it may be dominated by
        //few spikes
        findoutliers(upper,bins,&p01,&pos01, &p99,
                     &pos99, &binvalsmax, &binvalspos);
        if (pos99 <= 0)
        {
            binmax = binmin;
            delta=0.0;
            numbins = 0;
        }
        else
        {
            binmax = binmin + pos99 * bindelta;
            delta = width/pos99;
            numbins = pos99;
        }
    }
    else
    {
        if(binmax < -900)
        {
            findminmax_freq_only(upper,bins,&p00,&pos00, &p100,
                                 &pos100);
            if(pos100 <=0 )
            {
                binmax=binmin;
                delta=0.0;
                numbins=0;
            }
            else
            {
                binmax = binmin + pos100 * bindelta;
                delta = width/pos100;
                numbins = pos100;
            }
        }
    }
    if(valmax < binvalsmax)
    {
      binvalsmax = (double) valmax;
    }
    if(nofreq)
    {
        //	if(p99 > 5000.)
        //  p99 = 5000.;
        if(binvalsmax < p99)
            p99 = binvalsmax;
        
        pdf_create_axes(VERTICAL, 0, p99, 5, lx, ly, height);
    }
    else
    {
        pdf_create_axes(VERTICAL, 0, binvalsmax/total, 5, lx, ly, height);
    }
    pdf_create_axes(HORIZONTAL, binmin , binmax, 6, lx, ly, width);
    if(world->options->has_datefile)
    {
        ratio = width / (binmax - binmin);
        for(locus = 0; locus < world->loci; locus++)
        {
            if(world->bayes->mu)
                averagerate += world->options->meanmu[locus] * world->bayes->histogram[locus].modes[world->numpop2+world->locus];
            else
                averagerate += world->options->meanmu[locus];
        }
        averagerate /= world->loci;
        pop = topop;
        //	for(pop=0; pop < world->numpop; pop++)
        //  {
        if(world->data->numind != NULL)
        {
            for(ind = 0; ind < world->data->numind[pop][0]; ind++)
            {
	      xcoord = (double) (world->data->sampledates[pop][0][ind].date
				* world->options->generation_year
				* averagerate);
	      pdf_print_dot(lx+(xcoord * ratio), ly-3.0, 3.0, DIAMOND, blue);
            }
        }
	    //  }
    }
    
    
    //    x = delta/1.95;
    x = delta/2.0;
    
    xold = x;
    if(nofreq)
    {
      //nbmin = (long) MAX(2.,((double) numbins)/50.);//I need a better algorithm! MAX(2.,((double) numbins)/200.);
        //bayes_smooth(lower,numbins,nbmin, TRUE,FALSE);
        //bayes_smooth(binvals,numbins, nbmin, TRUE,FALSE);
        //bayes_smooth(upper,numbins, nbmin, TRUE,FALSE);
        warning("\nNO SMOOTHING for plotting!\n");
        vold = binvals[0];
    }
    for(i=0;i<numbins;i++)
    {
        if(set50[i] == '1')
            HPDF_Page_SetRGBStroke(page, 0.1f, 0.1f, 0.1f);
        else
        {
            if(set95[i] == '1')
                HPDF_Page_SetRGBStroke(page, 0.4f, 0.4f, 0.4f);
            else
                HPDF_Page_SetRGBStroke(page, 0.8f, 0.8f, 0.8f);
        }
        if(nofreq)
        {
            l = lower[i];
            u = upper[i];
            v = binvals[i];
            if(confidence !=NULL && confidence[i]>=1.0)
            {
                xold = x;
                vold = (double) v;
                x += delta;
                continue;
            }
            HPDF_Page_SetLineWidth(page, (float) delta);
            if(u/p99 > 1.0)
            {
                if(l/p99 < (1.0 + SMALLEPSILON))
                {
                    pdf_draw_line(lx+x,ly + ((double) l)/p99 * height,
                                  lx+x,ly + height);
                }
                pdf_print_dot(lx+x, ly + height, 3, SQUARE, red);
            }
            else
            {
                pdf_draw_line(lx+x,ly + ((double) l)/p99 * height,
                              lx+x,ly + ((double) u)/p99 * height);
            }
            HPDF_Page_SetRGBStroke(page, 0.2f, 0.2f, 0.2f);
            
            if(((double) v)/p99 < (1.0 + SMALLEPSILON))
            {
                //                pdf_draw_line(lx+x,ly + ((double) v)/p99 * height-1.0F,
                //              lx+x,ly + ((double) v)/p99 * height+1.0F);
                HPDF_Page_SetLineWidth(page, 1.0);
                pdf_draw_line(lx+xold,ly + ((double) vold)/p99 * height,
                              lx+x,ly + ((double) v)/p99 * height);
                vold  = (double) v;
                xold = x;
            }
            else
            {
                if(vold<p99)
                {
                    HPDF_Page_SetLineWidth(page, 1.0);
                    pdf_draw_line(lx+xold,ly + ((double) vold)/p99 * height,
                                  lx+x,ly + height);
                }
                vold  = p99;
                xold = x;
            }
        }
        else
        {
            HPDF_Page_SetRGBStroke(page, 0.8f, 0.8f, 0.8f);
            HPDF_Page_SetLineWidth(page, (float) delta);
            pdf_draw_line(lx+x,ly+0,lx+x,ly + ((double) binvals[i])/binvalsmax * height);
            HPDF_Page_SetLineWidth(page, (float) (delta/3.));
            HPDF_Page_SetRGBStroke(page, 0.5f, 0.5f, 0.5f);
            pdf_draw_line(lx+x,ly+((double) binvals[i])/binvalsmax * height,
                          lx+x,ly + ((double) upper[i])/binvalsmax * height);
            pdf_draw_line(lx+x,ly+((double) lower[i])/binvalsmax * height,
                          lx+x,ly + ((double) binvals[i])/binvalsmax * height);
        }
        //fprintf(stderr,"%f %f %f\n",lx+x,ly+ binvals[i]/binvalsmax * height, binvals[i]);
        if(confidence != NULL)
        {
	  r = (float) confidence[i];
	  if(r>1.0f && r < 0.0f)
	    continue;
	  HPDF_Page_SetLineWidth(page, (float) delta);
	  HPDF_Page_SetRGBStroke(page, r,r,r);
	  pdf_draw_line(lx+x,ly + height+5.0,
			lx+x,ly + height+10.0);
            //	    pdf_print_dot(lx+x, ly + height + 5, 3, SQUARE, red);
        }
        x += delta;
    }
    HPDF_Page_SetRGBStroke(page, 0.f, 0.f, 0.f);
    HPDF_Page_SetLineWidth(page, 1);
    myfree(lower);
    myfree(upper);
}

void pdf_master_init(world_fmt *world, option_fmt *options, data_fmt *data)
{
    ////////////////////////////////////
    pdf_init();
    // add first page
    pdf_new_page(options->title);
    pdf_master_title(options->title);//, &left_margin);
    pdf_print_options(world, options, data);
	if(options->datatype!='g')
    {
	    pdf_print_data_summary(world, options, data,  &page_height, &left_margin);
	    pdf_print_data (world, options, data);
	    if(options->murates_fromdata)
            pdf_print_mutationrate_weights(options->mu_rates, options->segregs, options->wattersons, world->loci);
    }
    ////////////////////////////////////
}


int pdf_init(void)
{
  //pdf_type1_fontdef font1_def;
  //pdf_type1_fontdef font2_def;
  //pdf_type1_fontdef font3_def;
  //pdf_type1_fontdef font4_def;
    
  /* Create a new PDF document. */
  //  doc = pdf_doc_new();
  //pdf_doc_new_doc(doc);
  doc = HPDF_New (error_handler, NULL);
  if (!doc) {
    printf ("error: cannot create PdfDoc object\n");
    return 1;
  }
  
  //  if (setjmp(env)) {
  //  HPDF_Free (doc);
  //  return 1;
  //}
  helv  = HPDF_GetFont (doc, "Helvetica", NULL);
  helvob  = HPDF_GetFont (doc, "Helvetica-Oblique", NULL);
  symb  = HPDF_GetFont (doc, "Symbol", NULL);
  cour  = HPDF_GetFont (doc, "Courier", NULL);
    /* Add Helvetica Font. */
  //   font1_def = pdf_create_type1_fontdef(PDF_FONT_HELVETICA);
  //pdf_doc_add_type1font(doc, font1_def, NULL, NULL);
  /* Add Helvetica-Oblique Font. */
  //font2_def = pdf_create_type1_fontdef(PDF_FONT_HELVETICA_OBLIQUE);
  //pdf_doc_add_type1font(doc, font2_def, NULL, NULL);
  /* Add Symbol Font. */
  //font3_def = pdf_create_type1_fontdef(PDF_FONT_SYMBOL);
  //pdf_doc_add_type1font(doc, font3_def, NULL, NULL);
  /* Add Courier Font. */
  //font4_def = pdf_create_type1_fontdef(PDF_FONT_COURIRE);
  //pdf_doc_add_type1font(doc, font4_def, NULL, NULL);
  return 0;
}

///
/// generate a new PDF page with a border rectangle and with page numbers and
/// title in top right corner and impressum at bottom left
int pdf_new_page(char *title)
{
    char stemp[LINESIZE];
    if(strlen(title)>1)
        strncpy(pdf_pagetitle,title,80);
    /* Add a page to the document. */
    page = HPDF_AddPage(doc);
    HPDF_Page_SetSize (page, HPDF_PAGE_SIZE_LETTER, HPDF_PAGE_PORTRAIT);
    page_counter += 1;
    /* Get the page object of the page. */
    //page = pdf_page_get_page(page);
    
    /*Set current font to "Helvetica" and set the size to 10. */
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
    
    HPDF_Page_SetLineWidth(page, 1);
    /* draw page rectangle 50pt from border*/
    HPDF_Page_Rectangle(page, 50, 50,
                           HPDF_Page_GetWidth(page) - 100,
                           HPDF_Page_GetHeight(page) - 110);
    HPDF_Page_Stroke(page);
    /* print the title of the analysis*/
    pdf_print_header(pdf_pagetitle);
    /* print the impressum at the bottome*/
    sprintf(stemp,"Migrate %s: (http://popgen.sc.fsu.edu) [program run on %s]",MIGRATEVERSION, pdf_time);
    HPDF_Page_SetFontAndSize (page, helv, 6.0);
    pdf_print_contents_at(50, 42, stemp);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
    return 0;
}



void pdf_print_mutationrate_weights(MYREAL *murates, long *segregs, MYREAL *wattersons, long loci)
{
    double offset[]={55.,150.,280.,380.};
    double page_width;
    //double page_height;
    left_margin=55.0;
    char title[]="Relative mutation rate among loci estimated from the data";
    long locus;
    //long row;
    //long rows = loci+2;
    //long cols = 4;
    double lx = left_margin;
    double ly;
    char st[100];
    double mumean=0.0;
    double segregmean=0.0;
    double wamean=0.0;
    pdf_print_section_title(&page_width, &page_height, title);
    ly  = page_height;
    pdf_print_line_element(lx, ly, offset[0], "Locus");
    pdf_print_line_element(lx, ly, offset[1]-10.0, "Relative");
    if(wattersons!=NULL)
    {
        pdf_print_line_element(lx, ly, offset[2]-7.0, "Watterson's");
        symbol_Theta(lx+offset[2],ly ,11, -1);
        pdf_print_line_element(lx, ly, offset[3]+10.0, "Segregating");
    }
    else
    {
        pdf_print_line_element(lx, ly, offset[2], "Number of alleles");
    }
    pdf_advance(&page_height);
    ly = page_height;
    pdf_print_line_element(lx, ly, offset[1], "mutation rate");
    if(wattersons!=NULL)
    {
        pdf_print_line_element(lx, ly, offset[2], "(per site)");
        pdf_print_line_element(lx, ly, offset[3], "sites");
    }
    pdf_advance(&page_height);
    //ly = page_height;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_advance(&page_height);
    for (locus=0; locus < loci; locus++)
    {
        HPDF_Page_SetFontAndSize (page, helv, 10.0);
        sprintf(st,"%5li",locus+1);
        pdf_print_line_element(lx, page_height, offset[0], st);
        pdf_print_line_element2(lx, page_height, offset[1], (double) murates[locus],5,5);
        mumean += (murates[locus] - mumean)/(locus+1);
        if(wattersons != NULL)
        {
	  pdf_print_line_element2(lx, page_height, offset[2], (double) wattersons[locus],2,8);
            sprintf(st,"%6li",segregs[locus]);
            pdf_print_line_element(lx, page_height, offset[3], st);
            wamean += (wattersons[locus] - wamean)/(locus+1);
            segregmean += (segregs[locus] - segregmean)/(locus+1);
        }
        else
        {
            sprintf(st,"%6li",segregs[locus]);
            pdf_print_line_element(lx, page_height, offset[2], st);
            segregmean += (segregs[locus] - segregmean)/(locus+1);
        }
        pdf_advance(&page_height);
    }
    sprintf(st,"%6s","All");
    pdf_print_line_element(lx, page_height, offset[0], st);
    pdf_print_line_element2(lx, page_height, offset[1], mumean,5,5);
    if (wattersons!=NULL)
    {
        pdf_print_line_element2(lx, page_height, offset[2], wamean,2,8);
        pdf_print_line_element2(lx, page_height, offset[3], segregmean,5,1);
    }
    else
        pdf_print_line_element2(lx, page_height, offset[2], segregmean,5,1);
}

///
/// start a new page and print the section title on top of page
void pdf_print_section_title(double *page_width, double * mypage_height, char *title)
{
    double w;
    // setup new page and title
    pdf_new_page("");
    HPDF_Page_SetFontAndSize (page, helvob, 18.0);
    w = (double) HPDF_Page_TextWidth(page, title);
    *mypage_height = (double) HPDF_Page_GetHeight(page) - 75.0;
    *page_width = (double) HPDF_Page_GetWidth(page);
    *mypage_height -= 20.0;
    pdf_print_contents_at((*page_width - w)/2.0, *mypage_height, title);
    HPDF_Page_SetRGBStroke(page, 0.f, 0.f, 0.f);
    *mypage_height -= 20.0;
    pdf_draw_line(50, *mypage_height, *page_width-50.0, *mypage_height);
    *mypage_height -= 20.0;
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
}

///
/// print title on every page with page number
int pdf_print_header(char *title)
{
    double w;
    double page_width;
    char *fulltitle;
    
    fulltitle = (char*) mycalloc(255,sizeof(char));
    /* Print the title of the page (with positioning center). */
    sprintf(fulltitle,"%s -- %i",title, page_counter);
    //printf("%s\n",fulltitle);
    w = (double) HPDF_Page_TextWidth(page, fulltitle);
    /* Start to print text. */
    HPDF_Page_BeginText(page);
    /* Move the position of the text to center */
    page_height = (double) HPDF_Page_GetHeight(page);
    page_width = (double) HPDF_Page_GetWidth(page);
    float xx =  (float) (page_width - w - 50);
    float yy = (float) ( page_height - 50);
    //HPDF_Page_MoveTo(page, xx, yy);
    /* Print title with pagenumber to rightadjusted */
    HPDF_Page_TextOut(page, xx, yy, fulltitle);
    /* Finish to print text. */
    HPDF_Page_EndText(page);
    
    myfree(fulltitle);
    
    return 0;
}

void pdf_print_contents_at(double x, double y, char *title)
{
    HPDF_Page_BeginText(page);
    HPDF_Page_TextOut(page,  (float) x, (float) y, title);
    HPDF_Page_EndText(page);
    //fprintf(stderr,"(%f, %f) %s\n",x,y,title);
}

///
/// print first title page
int pdf_master_title(char *title)
{
#ifdef MPI
    int cpus;
#endif
    double w;
    double page_width;
    char newtitle[LINESIZE];
    
    get_time (pdf_time, "%H:%M:%S");
    
    /* Move the position of the text to center */
    
    page_height = (double) HPDF_Page_GetHeight(page);
    page_width = (double) HPDF_Page_GetWidth(page);
    /* Print the title of the page (with positioning center). */
    HPDF_Page_SetFontAndSize (page, helvob, 24.0);
    
    MYSNPRINTF(newtitle,(size_t) (page_width-110)/10, "%s",title);
    w = (double) HPDF_Page_TextWidth(page, newtitle);
    //*left_margin = 55;
    /* Start to print text. */
    ////HPDF_Page_BeginText(page);
    /* Print title centered */
    page_height -= 100;
    pdf_print_contents_at((page_width - w)/2, page_height, newtitle);
    HPDF_Page_SetRGBStroke(page, 0.f, 0.f, 0.f);
    page_height -= 26;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    HPDF_Page_SetFontAndSize (page, helv, 12.0);
    page_height -= 24;
    pdf_print_contents_at(55, page_height, "POPULATION SIZE, MIGRATION, DIVERGENCE, ASSIGNMENT, HISTORY");
    pdf_advance(&page_height);
    pdf_print_contents_at(55, page_height, "Bayesian inference using the structured coalescent");
    pdf_advance(&page_height);
    pdf_printf(55, page_height, 'L', "Migrate-n version %s [%s]",MIGRATEVERSION, MIGRATESUBVERSION);
#ifdef BEAGLE
    pdf_advance(&page_height);
    pdf_printf(55, page_height, 'L',"  Using HMSBEAGLE likelihood calculators\n");
#endif
#ifdef AVX
    pdf_advance(&page_height);
    pdf_printf(55, page_height, 'L',"  Using Intel AVX (Advanced Vector Extensions)\n");
#endif

#ifdef MPI
    pdf_advance(&page_height);
    pdf_printf(55, page_height, 'L', "  Compiled for PARALLEL computer architectures\n");
    pdf_advance(&page_height);
    cpus = numcpu - 1;
    pdf_printf(55, page_height, 'L', "  One master and %i compute nodes are available.\n", cpus);
#endif
#ifdef THREAD
    pdf_advance(&page_height);
    pdf_printf(55, page_height, 'L',"  Compiled for a SYMMETRIC multiprocessors\n");
#endif
#ifdef GRANDCENTRAL
    pdf_advance(&page_height);
    pdf_printf(55, page_height, 'L',"  Compiled for a SYMMETRIC multiprocessors (Grandcentral)\n");
#endif
#ifdef FAST_EXP
    pdf_advance(&page_height);
    pdf_printf(55, page_height,  'L',"  Fast approximation to Exp() and Log() used\n");
#endif
    pdf_advance(&page_height);
    pdf_print_time(55, &page_height, "Program started at ");
    firstpage = page;
    timestampy = page_height;
    // here the end time stamp will be printed at the very end the run
    pdf_advance(&page_height);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    // print migrate logo
    HPDF_Page_SetLineWidth(page, 1);
    HPDF_Page_SetRGBStroke(page, 0.f, 0.f, 0.f);
    HPDF_Page_SetRGBFill(page, 0.99f, 0.f, 0.f);
    pdf_migrate_logo(-1400, -1435 - LETTERADJUST, 100.0);
    HPDF_Page_SetRGBFill(page, 0.f, 0.1f, 0.9f);
    pdf_migrate_logo(-1364, -1435  - LETTERADJUST, 100.0);
    HPDF_Page_SetLineWidth(page, 4);
    pdf_migrate_logo_lines(-1400, -1435 - LETTERADJUST, 100.0);
    page_height = page_height - 3*LINESTRETCH;
    return 0;
}

int pdf_write_file(option_fmt *options)
{
    HPDF_SaveToFile(doc, options->pdfoutfilename);
    return 0;
}

///
/// print elements of Bayesian table for character variables
double pdf_print_line_element(double lx, double ly, double offset, char *title)
{
    double w=0;
    w = (double) HPDF_Page_TextWidth(page, title);
    if(offset>=0)
        pdf_print_contents_at(lx+offset-w, ly, title);
    else
        pdf_print_contents_at(lx-offset, ly, title);
    return w;
}

///
/// print elements of Bayesian table for double variables
double pdf_print_line_element2(double lx, double ly, double offset, double value, int fmt1, int fmt2)
{
  char sign[]=" ";
    double w=0;
    char title[100];
    if (fabs(value)>(double) HUGE)
      {
	if (value<0)
	  sprintf(sign,"-");
	snprintf(title,sizeof(title),"%sHUGE",sign);
      }
    else
      snprintf(title,sizeof(title),"%*.*f",fmt1,fmt2,value);
    w = (double) HPDF_Page_TextWidth(page, title);
    if(offset>=0)
        pdf_print_contents_at(lx+offset-w, ly, title);
    else
        pdf_print_contents_at(lx-offset, ly, title);
    return w;
}


void pdf_print_line_theta(double lx, double ly, double offset, long j)
{
    char tempstring[100];
    char * thetatitle="Q";
    if(j < 0)
        sprintf(tempstring," ");
    else
        sprintf(tempstring,"%li",j+1);
    HPDF_Page_SetFontAndSize (page, symb, 11.0);
    pdf_print_contents_at(lx-offset,ly,thetatitle);
    HPDF_Page_SetFontAndSize (page, symb, 8.0);
    pdf_print_contents_at(lx-offset+10,ly-4,tempstring);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
}
void pdf_print_line_growth(double lx, double ly, double offset, long j)
{
    char tempstring[100];
    char * growthtitle="g";
    if(j < 0)
        sprintf(tempstring," ");
    else
        sprintf(tempstring,"%li",j+1);
    HPDF_Page_SetFontAndSize (page, helv, 11.);
    pdf_print_contents_at(lx-offset,ly,growthtitle);
    HPDF_Page_SetFontAndSize (page, symb, 8.0);
    pdf_print_contents_at(lx-offset+10,ly-4,tempstring);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
}

void pdf_print_line_mlf(double lx, double ly, double offset, long j)
{
    char tempstring[100];
    char * growthtitle="a";
    if(j < 0)
        sprintf(tempstring," ");
    else
        sprintf(tempstring,"%li",j+1);
    HPDF_Page_SetFontAndSize (page, symb, 11.0);
    pdf_print_contents_at(lx-offset,ly,growthtitle);
    HPDF_Page_SetFontAndSize (page, symb, 8.0);
    pdf_print_contents_at(lx-offset+10,ly-4,tempstring);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
}

void pdf_print_line_species(char *title, double lx, double ly, double offset, long frompop, long topop)
{
    char tempstring[100];
    char tostring[100];
    // title is either "D" or "s";
    if(topop < 0)
        sprintf(tostring,"+");
    else
        sprintf(tostring,"%li",topop+1);
    sprintf(tempstring,"%li->%s",frompop+1,tostring);
    HPDF_Page_SetFontAndSize (page, symb, 11.0);
    pdf_print_contents_at(lx-offset,ly,title);
    HPDF_Page_SetFontAndSize (page, symb, 8.0);
    pdf_print_contents_at(lx-offset+10,ly-4,tempstring);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
}

/// print rate line
void pdf_print_line_rate(double lx, double ly, double offset, long j, long exponent)
{
    char tempstring[100];
    char * rtitle="m";
    if(j < 0)
        sprintf(tempstring," ");
    else
        sprintf(tempstring,"%li",j+1);
    HPDF_Page_SetFontAndSize (page, symb, 10.0);
    pdf_print_contents_at(lx-offset,ly,rtitle); // print mu
    HPDF_Page_SetFontAndSize (page, symb, 8.0);
    pdf_print_contents_at(lx-offset+10,ly-4,tempstring);//print locus number as subscript
    
    HPDF_Page_SetFontAndSize (page, symb, 10.0);
    sprintf(tempstring,"[10");
    pdf_print_contents_at(lx-offset+13,ly,tempstring); //print the scale should look like this [x10-5]
    HPDF_Page_SetFontAndSize (page, symb, 8.0);
    sprintf(tempstring,"%li",exponent);
    pdf_print_contents_at(lx-offset+26,ly+4,tempstring); // print the exponent as superscript
    HPDF_Page_SetFontAndSize (page, symb, 10.0);
    sprintf(tempstring,"]");
    pdf_print_contents_at(lx-offset+35,ly,tempstring);
    
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
}

///
/// pretty print the M values
void  pdf_print_line_mig(char *migtitle, double lx, double mypage_height, double offset, long frompop, long topop)
{
    double w;
    char tostring[100];
    char tempstring[100];
    if(topop < 0)
        sprintf(tostring,"+");
    else
        sprintf(tostring,"%li",topop+1);
    sprintf(tempstring,"%li->%s",frompop+1,tostring);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
    pdf_print_contents_at(lx-offset, mypage_height,migtitle);
    w = (double) HPDF_Page_TextWidth(page, migtitle);
    HPDF_Page_SetFontAndSize (page, helv, 8.0);
    pdf_print_contents_at(lx-offset+w,mypage_height-4,tempstring);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
}


///
/// Print Bayesian table header
/// and report the width of the text in the headerline
void pdf_print_bayestableheader(double mypage_height, double myleft_margin, double right_margin, double *offset)
{
  (void) right_margin;
  double lx = myleft_margin;
  double ly = mypage_height;
    
    pdf_print_line_element(lx, ly, offset[0], "Locus");
    pdf_print_line_element(lx, ly, offset[1], "Parameter");
    pdf_print_line_element(lx, ly, offset[2], "2.5%");
    pdf_print_line_element(lx, ly, offset[3], "25.0%");
    pdf_print_line_element(lx, ly, offset[4], "Mode");
    pdf_print_line_element(lx, ly, offset[5], "75.0%");
    pdf_print_line_element(lx, ly, offset[6], "97.5%");
    pdf_print_line_element(lx, ly, offset[7], "Median");
    pdf_print_line_element(lx, ly, offset[8], "Mean");
}

// Parameter        2.5%%      25.0%%   median    75.0%%   97.5%%     mode     mean\n"

///
/// print Bayesian table
void pdf_print_bayestable(world_fmt *world)
{
    int fmt1=8;
    int fmt2=5;
    double mu;
    long lmu;
    double meanmu;
    long j0, j;
    long l;
    long size = get_numparam(world);
    long npp = size ;
    long nppspec = npp - world->grownum;
    double lx;
    double *offset;
    bayes_fmt * bayes = world->bayes;
    bayeshistogram_fmt * hist;
    species_fmt *s;
    
    page_height = (double) HPDF_Page_GetHeight(page);
    double page_width = (double) HPDF_Page_GetWidth(page);
    
    left_margin = 60;
    double right_margin = page_width - 60;
    
    //    double w;
    
    char st[100];
    char title[100] = "Bayesian Analysis: Posterior distribution table";
    long frompop;
    long topop;
    long locus;
    long end = world->loci > 1 ? world->loci + 1 : 1;
    long start = world->options->tersepdf ? end - 1 : 0;
    //column to to right-align the table columns
    offset = (double *) mycalloc(9,sizeof(double));
    offset[0] = -1; //left align
    offset[1] = -40;
    offset[2] = 135;
    offset[3] = 191;
    offset[4] = 247;
    offset[5] = 303;
    offset[6] = 359;
    offset[7] = 415;
    offset[8] = 471;
    pdf_print_section_title(&page_width, &page_height, title);
    //page_height -= 126;
    lx = left_margin;
    page_height -= 20;
    pdf_draw_line(left_margin, page_height, right_margin, page_height);
    pdf_advance(&page_height);
    pdf_print_bayestableheader(page_height, left_margin, right_margin, offset);
    pdf_advance(&page_height);
    pdf_draw_line(left_margin, page_height, right_margin, page_height);
    pdf_advance(&page_height);
    for(locus=start; locus < end; locus++)
    {
        if(world->data->skiploci[locus])
	  continue;
	else
	  {
            mu = 1.0; //used to adjust values for rate estimates for others this is 1.
            hist = &bayes->histogram[locus];
            if(locus == world->loci)
                strcpy(st,"  All ");
            else
                sprintf(st,"%5li ",locus + 1);
            
            for(j0=0; j0< size; j0++)
            {
	      if(shortcut(j0,world,&j))
		{
		  continue;
		}
	      pdf_print_line_element(lx, page_height, offset[0], st);
	      if(j < world->numpop)
                {
		  pdf_print_line_theta(lx, page_height, offset[1], j0);
		  fmt1 = 8;
		  fmt2 = 5;
                }
	      else if (j< world->numpop2)
		{
		  m2mm(j0,world->numpop,&frompop, &topop);
		  if(world->options->usem)
		    {
		      pdf_print_line_mig("M", lx, page_height, offset[1], frompop, topop);
		      fmt1 = 8;
		      if (strchr (SEQUENCETYPES, world->options->datatype))
			fmt2 = 1;
		      else
			fmt2 = 3;
		    }
		  else
		    {
		      pdf_print_line_mig("xNm", lx, page_height, offset[1], frompop, topop);
		      fmt1 = 8;
		      fmt2 = 5;
		    }
		}
	      else if (world->bayes->mu && j==world->numpop2)
		{
		  // rate modifier used
		  if(locus==world->loci && end>1)
		    {
		      meanmu = 0.;
		      for(l=0;l<world->loci;l++)
			{
			  meanmu += world->options->meanmu[l];
			}
		      meanmu /= world->loci;
		      lmu = (long) (floor( log10(hist->modes[j] * meanmu)));
		      mu = meanmu * pow(10. , -lmu);
		    }
		  else
		    {
		      lmu = (long) (floor( log10(hist->modes[j]*world->options->meanmu[locus])));
		      mu = world->options->meanmu[locus] * pow(10. , -lmu);
		    }
		  pdf_print_line_rate(lx, page_height, offset[1], -1, lmu);
		  fmt1 = 8;
		  fmt2 = 5;
		}
	      else if (world->has_speciation && j < nppspec)
		{
		  // speciation
		  s = get_which_species_model(j,world->species_model,world->species_model_size);
		  if (j == s->paramindex_mu)
		    {
		      pdf_print_line_mig("D", lx, page_height, offset[1], s->from, s->to);
		    }
		  else
		    {
		      pdf_print_line_mig("S", lx, page_height, offset[1], s->from, s->to);
		    }
		  fmt1 = 8;
		  fmt2 = 5;
		}		      
	      else if(world->has_growth && j < npp)
		{
		  long d=j0-npp;
		  while (d < world->options->growpops_numalloc && world->options->growpops[d]==0)
		    {
		      d++;
		    }
		  pdf_print_line_growth(lx, page_height, offset[1], d);
		  fmt1 = 8;
		  fmt2 = 5;			      
		}
	    
	      // reset the to helvetica from symbol or small helvetica
	      HPDF_Page_SetFontAndSize (page, helv, 10.0);
	      //if (j<world->numpop2 && world->options->custm2[j]!='d')
	      //  {
	      pdf_print_line_element2(lx, page_height, offset[2], (double) (mu * hist->cred95l[j]),fmt1,fmt2);
	      pdf_print_line_element2(lx, page_height, offset[3], (double) (mu * hist->cred50l[j]),fmt1,fmt2);
	      pdf_print_line_element2(lx, page_height, offset[4], (double) (mu * hist->modes[j]),fmt1,fmt2);
	      pdf_print_line_element2(lx, page_height, offset[5], (double) (mu * hist->cred50u[j]),fmt1,fmt2);
	      pdf_print_line_element2(lx, page_height, offset[6], (double) (mu * hist->cred95u[j]),fmt1,fmt2);
	      pdf_print_line_element2(lx, page_height, offset[7], (double) (mu * hist->medians[j]),fmt1,fmt2);
	      pdf_print_line_element2(lx, page_height, offset[8], (double) (mu * hist->means[j]),fmt1,fmt2);
	      page_height -= LINESTRETCH;
	      //  }
	      if(page_height < 55)
		{
		  pdf_new_page("");
		  page_height = (double) HPDF_Page_GetHeight(page);
		  page_width = (double) HPDF_Page_GetWidth(page);
		  page_height -= 50;
		  lx = left_margin;
		  page_height -= 20;
		  pdf_draw_line(left_margin, page_height, right_margin, page_height);
		  page_height -= 20;
		  pdf_print_bayestableheader(page_height, left_margin, right_margin, offset);
		  page_height -= 20;
		  pdf_draw_line(left_margin, page_height, right_margin, page_height);
		  pdf_advance(&page_height);
		}
	    }
	    pdf_draw_line(left_margin, page_height, right_margin, page_height);
	    pdf_advance(&page_height);
	  } // only when loci contains info
    } // over all loci
    myfree(offset);
    pdf_print_citation("Bayesian inference", world);
}


///
/// print out the acceptance ratios for all the different Bayesian updates
void
pdf_bayes_print_accept(world_fmt *world)
{
    char title[LINESIZE];
    double w;
    left_margin = 55;
    double page_width;
    long j0,j=0;             //used to loop over all parameters
    long topop    =0;   // indicator into the parameter vector, specifying originating population
    long frompop  =0;   // receiving population
    char *stempo;       // string variable holding print-string
    char *stemp;        // pointer to string, seems to be need to don't get MYREAL free warnings
    long trials   =0;   //
    long tc = get_numparam(world);//position of genealogy accept rates
    bayes_fmt *bayes = world->bayes;
    species_fmt *  s;
    stempo = (char *) mycalloc(LINESIZE,sizeof(char));
    stemp = stempo;
    
    sprintf(title,"Acceptance ratios for all parameters and the genealogies");
    pdf_new_page("");
    HPDF_Page_SetFontAndSize (page, helvob, 16.0);
    w = (double) HPDF_Page_TextWidth(page, title);
    
    /* Start to print text. */
    page_height = (double) HPDF_Page_GetHeight(page);
    page_width = (double) HPDF_Page_GetWidth(page);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    HPDF_Page_SetRGBStroke(page, 0.f, 0.f, 0.f);
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
    
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    
    // This needs more attention but will need more stuff to safe
    if(world->options->datatype == 'g')
    {
	    pdf_print_contents_at(left_margin, page_height, "not available with datatype=Genealogy");
	    pdf_advance(&page_height);
	    myfree(stempo);
	    return;
    }
    
    pdf_print_contents_at(left_margin, page_height, "Parameter");
    pdf_print_contents_at(250, page_height, "Accepted changes");
    pdf_print_contents_at(450, page_height, "Ratio");
    pdf_advance(&page_height);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_advance(&page_height);

    // population sizes
    for(j0=0; j0 < world->numpop; j0++)
    {
        if(!strchr("c", bayes->custm2[j0]))
        {
            j = world->bayes->map[j0][1];
            if (j>=0)
            {
                if((trials=world->trials_archive[j])>0)
                {
                    symbol_Theta(left_margin, page_height, 12, j0+1);
                    pdf_printf(250, page_height, 'L', "%8li/%-8li",world->accept_archive[j],trials);
                    pdf_printf(450, page_height, 'L', "%8.5f", (MYREAL) world->accept_archive[j]/trials);
                    pdf_advance(&page_height);
                }
            }
        }
    }
    // migration rates
    for(j0=world->numpop; j0 < world->numpop2; j0++)
    {
        if(!strchr("0c", bayes->custm2[j0]))
        {
            j = world->bayes->map[j0][1];
            if (j>=0)
            {
                if((trials=world->trials_archive[j])>0)
                {
                    m2mm (j0, world->numpop, &frompop, &topop);
                    symbol_M(left_margin, page_height, 12, frompop+1, topop+1, world->options->usem);
                    pdf_printf(250, page_height, 'L', "%8li/%-8li",world->accept_archive[j],trials);
                    pdf_printf(450, page_height, 'L', "%8.5f", (MYREAL) world->accept_archive[j]/trials);
                    pdf_advance(&page_height);
                    memset(stemp,0,sizeof(char)*(LINESIZE-1));
                }
            }
        }
    }
    // accepted rate of mutation rate changes
    if(bayes->mu)
    {
        if((trials=world->trials_archive[j0])>0)
        {
	  j=world->numpop2;
	  symbol_R(left_margin, page_height, 12, -1);
	  pdf_printf(250, page_height, 'L', "%8li/%-8li",world->accept_archive[j],trials);
	  pdf_printf(450, page_height, 'L', "%8.5f", (MYREAL) world->accept_archive[j]/trials);
	  pdf_advance(&page_height);
        }
    }
    // accepted speciation times and variances
    if(world->has_speciation)
    {
      for(j=0; j < world->species_model_size;j++)
        {
            s = &world->species_model[j];
            j0 = world->numpop2 + bayes->mu + 2 * s->id;
	    trials=world->trials_archive[j0];
            if(trials>0)
            {
                symbol_D(left_margin, page_height, 12, s->from+1,s->to+1);
                pdf_printf(250, page_height, 'L', "%8li/%-8li",world->accept_archive[j0],trials);
                pdf_printf(450, page_height, 'L', "%8.5f", (MYREAL) world->accept_archive[j0]/trials);
                pdf_advance(&page_height);
            }
	    trials=world->trials_archive[j0+1];
            if(trials>0)
            {
                symbol_S(left_margin, page_height, 12, s->from+1,s->to+1);
                pdf_printf(250, page_height, 'L', "%8li/%-8li",world->accept_archive[j0+1],trials);
                pdf_printf(450, page_height, 'L', "%8.5f", (MYREAL) world->accept_archive[j0+1]/trials);
                pdf_advance(&page_height);
            }
        }
    }
    if(world->has_growth)
      {
	for (j=0;j<world->grownum;j++)
	  {
	    j0 = world->numpop2 + bayes->mu + 2 * world->species_model_size + j;
	    trials=world->trials_archive[j0];
	    long d=j;
	    while (d < world->options->growpops_numalloc && world->options->growpops[d]==0)
	      d++;
	    symbol_Growth(left_margin, page_height, 12, d+1);
	    pdf_printf(250, page_height, 'L', "%8li/%-8li",world->accept_archive[j0],trials);
	    pdf_printf(450, page_height, 'L', "%8.5f", (MYREAL) world->accept_archive[j0]/trials);
	    pdf_advance(&page_height);
	    
	  }
      }
    // accepted trees
    if((trials=world->trials_archive[tc])>0)
    {
        pdf_printf(left_margin, page_height,'L', "Genealogies");
        pdf_printf(250, page_height, 'L', "%8li/%-8li", world->accept_archive[tc], trials);
        pdf_printf(450, page_height, 'L', "%8.5f", (MYREAL) world->accept_archive[tc]/ trials);
        pdf_advance(&page_height);
    }
    myfree(stempo);
}
///
/// print out the hyperprior information

void
pdf_bayes_print_hyperpriors(world_fmt *world)
{
    char title[LINESIZE];
    double w;
    left_margin = 55;
    double page_width;
    long j0,j=0;             //used to loop over all parameters
    long topop    =0;   // indicator into the parameter vector, specifying originating population
    long frompop  =0;   // receiving population
    char *stempo;       // string variable holding print-string
    char *stemp;        // pointer to string, seems to be need to don't get MYREAL free warnings
    long trials   =0;   //
    //long tc = world->numpop2 + world->bayes->mu + 2 * world->species_model_size; //position of genealogy accept rates
    bayes_fmt *bayes = world->bayes;
    hyper_fmt *hyperp = bayes->hyperp;
    species_fmt *  s;
    if (!bayes->hyperprior)
      return;
    stempo = (char *) mycalloc(LINESIZE,sizeof(char));
    stemp = stempo;
    sprintf(title,"Hyperpriors");
    pdf_new_page("");
    HPDF_Page_SetFontAndSize (page, helvob, 16.0);
    w = (double) HPDF_Page_TextWidth(page, title);
    
    /* Start to print text. */
    page_height = (double) HPDF_Page_GetHeight(page);
    page_width = (double) HPDF_Page_GetWidth(page);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    HPDF_Page_SetRGBStroke(page, 0.f, 0.f, 0.f);
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
    
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    
    // This needs more attention but will need more stuff to safe
    if(world->options->datatype == 'g')
    {
	    pdf_print_contents_at(left_margin, page_height, "not available with datatype=Genealogy");
	    pdf_advance(&page_height);
	    myfree(stempo);
	    return;
    }
    
    pdf_print_contents_at(left_margin, page_height, "Parameter");
    pdf_print_contents_at(250, page_height, "Priormean/std");
    pdf_print_contents_at(400, page_height, "PriorAlpha/std");
    pdf_print_contents_at(530, page_height, "N");
    pdf_advance(&page_height);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_advance(&page_height);
    // population sizes
    for(j0=0; j0 < world->numpop; j0++)
    {
      j = world->bayes->map[j0][1];
      if((trials=hyperp[j].meann)>1)
	{
	  symbol_Theta(left_margin, page_height, 12, j+1);
	  pdf_printf(250, page_height, 'L', "%8.5f/%-8.5f",hyperp[j].mean,hyperp[j].meanstd);
	  pdf_printf(400, page_height, 'L', "%8.5f/%-8.5f", hyperp[j].alpha,hyperp[j].alphastd);
	  pdf_printf(500, page_height, 'L', "%8li",trials);
	  pdf_advance(&page_height);
	}
    }
    // migration rates
    for(j0=world->numpop; j0 < world->numpop2; j0++)
    {
      j = world->bayes->map[j0][1];
      if (j>=0)
	{
	  trials=hyperp[j].meann;
	  if(trials>1)
	    {
	      m2mm (j0, world->numpop, &frompop, &topop);
	      symbol_M(left_margin, page_height, 12, frompop+1, topop+1, world->options->usem);
	      pdf_printf(250, page_height, 'L', "%8.5f/%-8.5f",hyperp[j].mean,hyperp[j].meanstd);
	      pdf_printf(400, page_height, 'L', "%8.5f/%-8.5f", hyperp[j].alpha,hyperp[j].alphastd);
	      pdf_printf(500, page_height, 'L', "%8li",trials);
	      pdf_advance(&page_height);
	      memset(stemp,0,sizeof(char)*(LINESIZE-1));
	    }
	}
    }
    // accepted rate of mutation rate changes
    if(bayes->mu)
    {
      j = world->bayes->map[world->numpop2][1];
      trials=hyperp[j].meann;
      if(trials>1)
        {
	  j=world->numpop2;
	  symbol_R(left_margin, page_height, 12, -1);
	  pdf_printf(250, page_height, 'L', "%8.5f/%-8.5f",hyperp[j].mean,hyperp[j].meanstd);
	  pdf_printf(400, page_height, 'L', "%8.5f/%-8.5f", hyperp[j].alpha,hyperp[j].alphastd);
	  pdf_printf(500, page_height, 'L', "%8li",trials);
	  pdf_advance(&page_height);
        }
    }
    // accepted speciation times and variances
    if(world->has_speciation)
    {
      for(j=0; j < world->species_model_size;j++)
        {
            s = &world->species_model[j];
            j0 = world->numpop2 + bayes->mu + 2 * s->id;
	    trials=hyperp[j0].meann;
            if(trials>1)
	      {
                symbol_D(left_margin, page_height, 12, s->from+1,s->to+1);
		pdf_printf(250, page_height, 'L', "%8.5f/%-8.5f",hyperp[j0].mean,hyperp[j0].meanstd);
		pdf_printf(400, page_height, 'L', "%8.5f/%-8.5f", hyperp[j0].alpha,hyperp[j0].alphastd);
		pdf_printf(500, page_height, 'L', "%8li",trials);
		pdf_advance(&page_height);
	      }

	    trials=hyperp[j0+1].meann;
            if(trials>1)
	      {
                symbol_S(left_margin, page_height, 12, s->from+1,s->to+1);
		pdf_printf(250, page_height, 'L', "%8.5f/%-8.5f",hyperp[j0+1].mean,hyperp[j0+1].meanstd);
		pdf_printf(400, page_height, 'L', "%8.5f/%-8.5f", hyperp[j0+1].alpha,hyperp[j0+1].alphastd);
		pdf_printf(500, page_height, 'L', "%8li",trials);
		pdf_advance(&page_height);
	      }
        }
    }
    myfree(stempo);
}

///
/// print out the autocorrelation in replicates and the effective sample size
void
pdf_bayes_print_ess(world_fmt *world)
{
    char title[LINESIZE];
    double w;
    left_margin = 55;
    double page_width;
    long j0=0;
    long j1=0;
    long j=0;             //used to loop over all parameters
    long topop    =0;   // indicator into the parameter vector, specifying originating population
    long frompop  =0;   // receiving population
    char *stempo;       // string variable holding print-string
    char *stemp;        // pointer to string, seems to be need to don't get MYREAL free warnings
    bayes_fmt *bayes = world->bayes;
    stempo = (char *) mycalloc(LINESIZE,sizeof(char));
    stemp = stempo;
    long tc = world->numpop2+world->bayes->mu + 2 * world->species_model_size;
    long tcg = world->numpop2+world->bayes->mu + 2 * world->species_model_size + world->grownum;
    sprintf(title,"MCMC-Autocorrelation and Effective MCMC Sample Size");
    pdf_new_page("");
    HPDF_Page_SetFontAndSize (page, helvob, 16.0);
    w = (double) HPDF_Page_TextWidth(page, title);
    
    /* Start to print text. */
    page_height = (double) HPDF_Page_GetHeight(page);
    page_width = (double) HPDF_Page_GetWidth(page);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    HPDF_Page_SetRGBStroke(page, 0.f, 0.f, 0.f);
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
    
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    
    pdf_print_contents_at(left_margin, page_height, "Parameter");
    pdf_print_contents_at(250, page_height, "Autocorrelation");
    pdf_print_contents_at(450, page_height, "Effective Sampe Size");
    pdf_advance(&page_height);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_advance(&page_height);
    // population sizes
    for(j0=0; j0 < world->numpop; j0++)
      {
	if(!strchr("c", world->options->custm2[j0]))
	  {
	    j = world->bayes->map[j0][1];
	    symbol_Theta(left_margin, page_height, 12, j0+1);
	    pdf_printf(250, page_height, 'L', "%5.5f",world->auto_archive[j]);
	    pdf_printf(450, page_height, 'L', "%10.2f", world->ess_archive[j]);
	    pdf_advance(&page_height);
	  }
      }
    // migration rates
    for(j0=world->numpop; j0 < world->numpop2; j0++)
      {
	if(!strchr("0c", bayes->custm2[j0]))
	  {
	    j = world->bayes->map[j0][1];
	    if (j>=0)
	      {
		m2mm (j0, world->numpop, &frompop, &topop);
		symbol_M(left_margin, page_height, 12, frompop+1, topop+1, world->options->usem);
		pdf_printf(250, page_height, 'L', "%5.5f",world->auto_archive[j]);
		pdf_printf(450, page_height, 'L', "%10.2f", world->ess_archive[j]);
		pdf_advance(&page_height);
	      }
	    memset(stemp,0,sizeof(char)*(LINESIZE-1));
	  }
      }
    j = world->numpop2; //adjusting j to fit
    // accepted rate of mutation rate changes
    if(world->bayes->mu)
    {
      symbol_R(left_margin, page_height, 12, -1);
      pdf_printf(250, page_height, 'L', "%5.5f",world->auto_archive[j]);
      pdf_printf(450, page_height, 'L', "%10.2f", world->ess_archive[j]);
      pdf_advance(&page_height);
    }
    if(world->species_model_size > 0)
      {
	for(j0=world->numpop2+world->bayes->mu;
	    j0<tc;
	    j0++)
	  {
	    if(shortcut(j0,world,&j))
	      {
		continue;
	      }
	    else
	      {
		species_fmt *s = get_which_species_model(j,world->species_model,world->species_model_size);
		frompop = s->from;
		topop = s->to;
		if(j == s->paramindex_mu)
		  {
		    symbol_D(left_margin, page_height, 12, frompop+1, topop+1);
		    pdf_printf(250, page_height, 'L', "%5.5f",world->auto_archive[j]);
		    pdf_printf(450, page_height, 'L', "%10.2f", world->ess_archive[j]);
		    pdf_advance(&page_height);
		  }
		else
		  {
		    symbol_S(left_margin, page_height, 12, frompop+1, topop+1);
		    pdf_printf(250, page_height, 'L', "%5.5f",world->auto_archive[j]);
		    pdf_printf(450, page_height, 'L', "%10.2f", world->ess_archive[j]);
		    pdf_advance(&page_height);		    
		  }
	      }
	  }
      }
    // accepted growth
    if(world->has_growth)
    {
        for (j1=tc; j1<tc+world->grownum; j1++)
        {
	  long d=j1-tc;
	  while (d < world->options->growpops_numalloc && world->options->growpops[d]==0)
	    d++;
	  symbol_Growth(left_margin, page_height, 12, d+1);
	  pdf_printf(250, page_height, 'L', "%5.5f",world->auto_archive[j1]);
	  pdf_printf(450, page_height, 'L', "%10.2f", world->ess_archive[j1]);
	  pdf_advance(&page_height);
        }
    }
    // accepted growth
    // likelihood of trees
    pdf_print_contents_at(left_margin, page_height,"Genealogies");
    pdf_printf(250, page_height, 'L', "%5.5f",world->auto_archive[j1]);
    pdf_printf(450, page_height, 'L', "%10.2f", world->ess_archive[j1]);
    pdf_advance(&page_height);
    
    myfree(stempo);
}


///
/// print out marginal likelihood calculated from thermodynamic integration, steppingstone sampling, and harmonic mean
void pdf_bayes_factor_header(world_fmt *world, option_fmt *options)
{
  (void) world;
  (void) options;
    char title[LINESIZE];
    double w;
    left_margin = 55;
    double page_width;
    sprintf(title,"Log-Probability of the data given the model (marginal likelihood)");
    pdf_new_page("");
    HPDF_Page_SetFontAndSize (page, helvob, 16.0);
    w = (double) HPDF_Page_TextWidth(page, title);
    
    /* Start to print text. */
    page_height = (double) HPDF_Page_GetHeight(page);
    page_width = (double) HPDF_Page_GetWidth(page);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    HPDF_Page_SetRGBStroke(page, 0.f, 0.f, 0.f);
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
    
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height, 'L', "Use this value for Bayes factor calculations:");
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height,'L', "BF = Exp[ ln(Prob(D | thisModel) - ln( Prob( D | otherModel)");
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height,'L', "or as LBF = 2 (ln(Prob(D | thisModel) - ln( Prob( D | otherModel))");
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height,'L', "shows the support for thisModel]");
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_advance(&page_height);
}


void pdf_bayes_factor_rawscores_header(world_fmt *world, option_fmt *options)
{
  (void) world;
  (void) options;
    left_margin = 55;
    double page_width;
    page_width = (double) HPDF_Page_GetWidth(page);
    pdf_advance(&page_height);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_advance(&page_height);
    pdf_print_contents_at(left_margin, page_height, "Locus");
    pdf_print_contents_at(170, page_height, "TI(1a)");
    pdf_print_contents_at(270, page_height, "BTI(1b)");
    //pdf_print_contents_at(380, page_height, "SS(2)");
    //pdf_print_contents_at(495, page_height, "HS(3)");
    pdf_print_contents_at(380, page_height, "HS(3)");
    pdf_advance(&page_height);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_advance(&page_height);
}

void pdf_bayes_factor_rawscores(long locus, MYREAL rawtermo, MYREAL beziertermo, MYREAL ss, MYREAL harmo)
{
  (void) ss;

    double page_width;
    page_width = (double) HPDF_Page_GetWidth(page);
    if(locus<0)
    {
        pdf_draw_line(50, page_height,  page_width-50, page_height);
        pdf_advance(&page_height);
        pdf_printf_right(left_margin+20, page_height,"All  ");
    }
    else
        pdf_printf_right(left_margin+20, page_height,"%5li", locus+1);
    pdf_printf_right(200, page_height,"   %12.2f", rawtermo);
    pdf_printf_right(303, page_height,"   %12.2f", beziertermo);
    //pdf_printf_right(407, page_height,"   %12.2f", ss);
    //pdf_printf_right(520, page_height,"   %12.2f", harmo);
    pdf_printf_right(407, page_height,"   %12.2f", harmo);
    pdf_advance(&page_height);
}
void pdf_bayes_factor_rawscores_minmax(int minmax, MYREAL rawtermo, MYREAL beziertermo, double harmo)
{
    double page_width;
    page_width = (double) HPDF_Page_GetWidth(page);
    pdf_draw_line(50, page_height,  page_width-50, page_height);
    pdf_advance(&page_height);
    if (minmax<0)
      pdf_printf_right(left_margin+20, page_height,"Lowest");
    else
      pdf_printf_right(left_margin+20, page_height,"Highest");
    pdf_printf_right(200, page_height,"   %12.2f", rawtermo);
    pdf_printf_right(303, page_height,"   %12.2f", beziertermo);
    //pdf_printf_right(407, page_height,"   %12.2f", ss);
    //pdf_printf_right(520, page_height,"   %12.2f", harmo);
    pdf_printf_right(407, page_height,"   %12.2f", harmo);
    pdf_advance(&page_height);
}

/*void pdf_bayes_factor_rawscores_harmo(long locus, MYREAL harmo)
{
    double page_width;
    page_width = HPDF_Page_GetWidth(page);
    if(locus < 0)
    {
        pdf_draw_line(50, page_height, page_width-50, page_height);
        pdf_advance(&page_height);
        pdf_printf_right(left_margin+20, page_height,"All  ");
    }
    else
        pdf_printf(left_margin, page_height,'R',"%5li", locus+1);
    pdf_draw_line(140, page_height, 150, page_height);
    pdf_draw_line(280, page_height, 290, page_height);
    pdf_printf_right(520, page_height,"   %20.2f", harmo);
    pdf_advance(&page_height);
    }*/

///
/// print out marginal likelihood calculated from thermodynamic integration, harmonic and arithmetic mean.
void
pdf_bayes_factor_comment(world_fmt *world,  MYREAL scaling_factor)
{
    left_margin = 55;
    //double page_width;
    //page_width = HPDF_Page_GetWidth(page);
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height,'L',"(1a) TI: Thermodynamic integration: log(Prob(D|Model)): Good approximation with many temperatures\n");
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height,'L',"(1b) BTI: Bezier-approximated Thermodynamic integration: when using few temperatures USE THIS!\n");
    //pdf_advance(&page_height);
    //pdf_printf(left_margin, page_height,'L',"(2)  SS: Steppingstone Sampling (Xie et al 2011)\n");
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height,'L',"(3)  HS: Harmonic mean approximation: Overestimates the marginal likelihood, poor variance\n\n");
    pdf_advance(&page_height);
    switch(world->options->adaptiveheat)
    {
        case STANDARD:
            pdf_printf(left_margin, page_height,'L',"%s","Adaptive heating was ON, therefore the values of (1) may be incorrect),");
            pdf_advance(&page_height);
            break;
        case BOUNDED:
            pdf_printf(left_margin, page_height,'L',"%s","Adaptive heating with bounds was ON, therefore the values of (1) may be incorrect),");
            pdf_advance(&page_height);
            break;
    }
    if(world->loci>1)
    {
        pdf_printf(left_margin, page_height,'L',"[Scaling factor = %f]", scaling_factor);
    }
    pdf_advance(&page_height);
    pdf_print_citation( "Marginal likelihood", world);
}

void pdf_burnin_stops(world_fmt *world, long maxreplicate)
{
    long z;
    double w;
    char title[LINESIZE];
    left_margin = 55;
    double page_width;
    sprintf(title,"Stop of burnin-in phase due to convergence");
    pdf_new_page("");
    HPDF_Page_SetFontAndSize (page, helvob, 16.0);
    w = (double) HPDF_Page_TextWidth(page, title);
    
    /* Start to print text. */
    page_height = (double) HPDF_Page_GetHeight(page);
    page_width = (double) HPDF_Page_GetWidth(page);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    HPDF_Page_SetRGBStroke(page, 0.f, 0.f, 0.f);
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
    
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    pdf_print_contents_at(left_margin, page_height, "Locus");
    pdf_print_contents_at(120, page_height, "Replicate");
    pdf_print_contents_at(210, page_height, "Steps");
    pdf_print_contents_at(310, page_height, "Variance ratio (new/old variance)");
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_advance(&page_height);
    for(z=0; z < world->loci * maxreplicate; z++)
    {
        pdf_printf(left_margin, page_height,'L',"%5li", world->burnin_stops[z].locus);
        pdf_printf(120, page_height,'L',"%10li",world->burnin_stops[z].replicate);
        pdf_printf(200, page_height,'L',"%10li", world->burnin_stops[z].stopstep);
        pdf_printf(310, page_height,'L',"%f", world->burnin_stops[z].variance/world->burnin_stops[z].oldvariance);
        pdf_printf(390, page_height,'L',"(%f/%f)",world->burnin_stops[z].variance,
                   world->burnin_stops[z].oldvariance);
        pdf_advance(&page_height);
    }
    pdf_advance(&page_height);
    pdf_advance(&page_height);
}

void pdf_print_stored_warnings(world_fmt *world)
{
    double w;
    char title[LINESIZE];
    left_margin = 55;
    double page_width;
    char *buffer;
    char *b;
    char *tmp;
    long i=0;
    long z=0;
    char *section;
    char paragraph[] = "This section reports potential problems with your run, but such \
reporting is often not very accurate. With many parameters in a \
multilocus analysis, it is common that some parameters for some \
loci will not be informative. These parameters then trigger suggestions \
to increase the prior range that are not sensible. Do not blindly \
follow the suggestions given. If some parameters are flagged, inspect \
the tables carefully and judge whether an action is required. Suppose \
you run a Bayesian inference with sequence data for macroscopic \
species. In that case, there is rarely the need to increase the \
prior for Theta beyond 0.1. If you use microsatellites data, it is \
rather common that your prior distribution for Theta should have a \
range from 0.0 to 100 or more. With many populations (>3), it is \
also very common that some migration routes are estimated poorly \
because the data contains little or no information for that route. \
Increasing the prior range will not help in such situations, but \
reducing the number of parameters may help.";
    section = (char *) mycalloc(LINESIZE,sizeof(char));
    sprintf(title,"Potential Problems");
    pdf_new_page("");
    HPDF_Page_SetFontAndSize (page, helvob, 16.0);
    w = (double) HPDF_Page_TextWidth(page, title);
    
    /* Start to print text. */
    page_height = (double) HPDF_Page_GetHeight(page);
    page_width = (double) HPDF_Page_GetWidth(page);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    HPDF_Page_SetRGBStroke(page, 0.f, 0.f, 0.f);
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
    
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    //w = 0;
    while(paragraph[i]!='\0')
    {
        section[z]=paragraph[i];
        section[z+1] = '\0';
        z++;
        i++;
        w = (double) HPDF_Page_TextWidth(page, section);
        if(w >= page_width-110)
        {
            while(section[z-1]!=' ')
            {
                z--;
                i--;
            }
            section[z]='\0';
            pdf_printf_next(left_margin, &page_height,section);
            section[0]='\0';
            z=0;
        }
    }
    pdf_printf_next(left_margin, &page_height,section);
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    if(world->warningsize > 0)
    {
        buffer = (char *) mycalloc(strlen(world->warning)+1,sizeof(char));
        sprintf(buffer,"%s",world->warning);
        b = buffer;
        tmp = strsep(&buffer,"\n");
        while(tmp!=NULL)
        {
            pdf_printf(left_margin, page_height,'L',"%s",tmp);
            pdf_advance(&page_height);
            tmp = strsep(&buffer,"\n");
        }
        myfree(b);
    }
    else
    {
        pdf_printf(left_margin, page_height,'L',"No warning was recorded during the run");
    }
    myfree(section);
}


///
/// fill and stroke bezier curve
void pdf_fill_stroke(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
  HPDF_Page_MoveTo(page, (float) x1, (float) y1);
    HPDF_Page_CurveTo(page, (float) x2, (float) y2, (float) x3, (float) y3, (float) x4, (float) y4);
    HPDF_Page_FillStroke(page);
}

///
/// plot migrate logo
void pdf_migrate_logo(double x, double y, double stretch)
{
  HPDF_Page_MoveTo(page, (float) (x + stretch * 18.705982), (float) (y + stretch * 21.16717));
  HPDF_Page_CurveTo(page, (float)  (x + stretch * 18.705982),
		    (float) (y + stretch * 21.234596), (float) (x + stretch * 18.779981),(float) (y + stretch * 21.289254),
		    (float) (x + stretch * 18.871266), (float) (y + stretch * 21.289254));
  HPDF_Page_CurveTo(page, (float) (x + stretch * 18.962551),
		    (float) (y + stretch * 21.289254), (float) (x + stretch * 19.03655), (float) (y + stretch * 21.234596),
		    (float) (x + stretch * 19.03655), (float) (y + stretch * 21.167171));
  HPDF_Page_CurveTo(page, (float) (x + stretch * 19.03655),
		    (float) (y + stretch *  21.099745), (float) (x + stretch * 18.962551), (float) (y + stretch * 21.045087),
		    (float) (x + stretch * 18.871266), (float) (y + stretch * 21.045087));
  HPDF_Page_CurveTo(page, (float) (x + stretch * 18.779981),
		    (float) (y + stretch * 21.045087), (float) (x + stretch * 18.705982), (float) (y + stretch * 21.099745),
		    (float) (x + stretch * 18.705982), (float) (y + stretch * 21.167171));
    HPDF_Page_FillStroke(page);
}

void pdf_migrate_logo_lines(double x, double y, double stretch)
{
  HPDF_Page_MoveTo(page, (float) (x + stretch * 18.773599), (float) (y + stretch * 21.177612));
  HPDF_Page_LineTo(page, (float) (x + stretch * 18.773599), (float) (y + stretch * 20.997156));
  HPDF_Page_LineTo(page, (float) (x + stretch * 18.882535), (float) (y + stretch * 20.997156));
  HPDF_Page_LineTo(page, (float) (x + stretch * 18.882861), (float) (y + stretch * 21.177612));
  HPDF_Page_Stroke(page);
  HPDF_Page_MoveTo(page, (float) (x + stretch * 18.979539),  (float) (y + stretch * 21.177612));
  HPDF_Page_LineTo(page, (float) (x + stretch * 18.980202),  (float) (y + stretch * 20.948318));
    HPDF_Page_LineTo(page, (float) (x + stretch * 19.104166),  (float) (y + stretch * 20.948318));
    HPDF_Page_LineTo(page, (float) (x + stretch * 19.10341 ),  (float) (y + stretch * 21.177612));
    HPDF_Page_Stroke(page);
    HPDF_Page_MoveTo(page, (float) (x + stretch * 19.213683),  (float) (y + stretch * 21.177612));
    HPDF_Page_LineTo(page, (float) (x + stretch * 19.213103),  (float) (y + stretch * 20.891974));
    HPDF_Page_LineTo(page, (float) (x + stretch * 19.045865),  (float) (y + stretch * 20.891974));
    HPDF_Page_LineTo(page, (float) (x + stretch * 19.045865),  (float) (y + stretch * 20.948318));
    HPDF_Page_Stroke(page);					
    HPDF_Page_MoveTo(page, (float) (x + stretch * 19.318285),  (float) ( y + stretch * 21.177612));
    HPDF_Page_LineTo(page, (float) (x + stretch * 19.318285),  (float) ( y + stretch * 20.809329));
    HPDF_Page_LineTo(page, (float) (x + stretch * 19.132266),  (float) ( y + stretch * 20.809329));
    HPDF_Page_LineTo(page, (float) (x + stretch * 19.132521),  (float) ( y + stretch * 20.890561));
    HPDF_Page_Stroke(page);					
    HPDF_Page_MoveTo(page, (float) (x + stretch * 19.228543),  (float) (y + stretch * 20.645554));
    HPDF_Page_LineTo(page, (float) (x + stretch * 19.229199),  (float) (y + stretch * 20.808985));
    HPDF_Page_Stroke(page);					
    HPDF_Page_MoveTo(page, (float) (x + stretch * 18.829904),  (float) (y + stretch * 20.647076));
    HPDF_Page_LineTo(page, (float) (x + stretch * 18.829904),  (float) (y + stretch * 20.996422));
    HPDF_Page_Stroke(page);
    
}

///
/// print time stamp in pdf
void pdf_print_time(double startx, double *pageheight, char text[])
{
    char nowstr[LINESIZE];
    get_time(nowstr, "  %c");
    if (nowstr[0] != '\0')
        pdf_printf(startx, *pageheight, 'L',"%s %s", text, nowstr);
    pdf_advance(pageheight);
}

///
/// prints time stamp for end of run needs pretty.c-global variables firstpage and timestampy/// to successfully print timestamp, compare to call of pdf_print_time() function.
void pdf_print_end_time(double *mypage_height)
{
  (void) mypage_height;
    HPDF_Page thispage = firstpage;
    char nowstr[LINESIZE];
    char *title;
    time_t endseconds;
    char runtime[LINESIZE];
    title = (char *) mycalloc(LINESIZE,sizeof(char));
    
    get_time(nowstr, "  %c");
    endseconds = time(0);
    get_runtime(runtime,startseconds,endseconds);
    sprintf(title,"Program finished at %s [%s]", nowstr, runtime);
    
    if (nowstr[0] != '\0')
    {
        HPDF_Page_SetFontAndSize (thispage, helv, 12.0);
        HPDF_Page_BeginText(thispage);
	HPDF_Page_TextOut(thispage, 55., (float) timestampy, title);
        HPDF_Page_EndText(thispage);
    }
    myfree(title);
}

///
/// add a new page and set the master title
void pdf_title(char *title, double page_width)
{
  double ph = page_height;
    double w;
    pdf_new_page("");
    HPDF_Page_SetFontAndSize (page, helvob, 16.0);
    w = (double) HPDF_Page_TextWidth(page, title);
    pdf_print_contents_at((page_width - w)/2, ph - 100, title);
    HPDF_Page_SetRGBStroke(page, 0.f, 0.f, 0.f);
    pdf_draw_line(50, ph - 126, page_width-50, ph - 126);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
    page_height -= 126.0;
}

void find_posterior_min_max(double *minval, double *maxval, long start, long stop, bayes_fmt *bayes, long locus)
{
    long j, pa0,pa;
    long numbins  = 0;
    long *bins = bayes->histogram[locus].bins;
    double *results = bayes->histogram[locus].results;
    MYREAL *mini = bayes->histogram[locus].minima;
    MYREAL *maxi = bayes->histogram[locus].maxima;
    double p01, p99, mm;
    double val;
    long pos01, pos99, posmm;
    // set all parameter print maxima to the same value for theta and M
    *minval  = (double) HUGE;
    *maxval  = 0.;
    for(pa=0;pa < start; pa++)
        numbins += bins[pa];
    for(pa0=start; pa0 < stop; pa0++)
    {
        // if custom migration matrix is set to zero
        // continue
        //if(strchr("0c",bayes->custm2[pa]))
        //continue;
        if(bayes->map[pa0][1] == INVALID)
	  continue;
        else
            pa = bayes->map[pa0][1];
        for(j=0;j < pa; j++)
            numbins += bins[j];
        findmoments(&results[numbins],bins[pa],&p01,&pos01,&p99,&pos99,&mm,&posmm);
        if(mini[pa] < *minval)
	  *minval = (double) mini[pa];
        val = (double) (pos99 * (maxi[pa]-mini[pa])/bins[pa]);
        if(val> *maxval)
            *maxval = val;
        //      numbins += bins[pa];
    }
}


void pdf_pretty_histogram(long pa, long rpa,long numbins, double * results, long stride, char * set50,
                          char * set95, long * bins, double delta, MYREAL * mini, MYREAL * maxi,
                          double lx, double ly, world_fmt *world, double themin, double themax)
{
  (void) stride;
    bayes_fmt *bayes = world->bayes;
    long p99bins;
    double minival;
    double maxival;
    //printf("%li -------@@\n",pa);
    if (mini==NULL)
        minival = themin;
    else
      minival = (double) mini[rpa];
    if (maxi==NULL)
        maxival = themax;
    else
      maxival = (double) maxi[rpa];
    switch(bayes->prettyhist)
    {
        case PRETTY_MAX:
            pdf_histogram( &results[numbins],&set50[numbins], &set95[numbins],
                          bins[pa],(double) delta, (double) minival, (double) maxival,
                          lx,ly,187,116, FALSE, &world->bayes->priors[numbins]);
            break;
        case PRETTY_P99:
            pdf_histogram( &results[numbins],&set50[numbins], &set95[numbins],
                          bins[rpa],(double) delta, (double) minival, -9999,
                          lx,ly,187,116, FALSE, &world->bayes->priors[numbins]);
            break;
        case PRETTY_P99MAX:
            p99bins = (long)((themax-themin)/(double)delta);
            pdf_histogram( &results[numbins],&set50[numbins], &set95[numbins],
                          p99bins,(double) delta, themin, themax,
                          lx,ly,187,116, FALSE, &world->bayes->priors[numbins]);
            break;
        case PRETTY_P100:
        default:
            pdf_histogram( &results[numbins],&set50[numbins], &set95[numbins],
                          bins[rpa],(double) delta, (double) minival, -999,
                          lx,ly,187,116, FALSE, &world->bayes->priors[numbins]);
            break;
    }
}


void pdf_advance_hist(long numpop, long *z, double *lx, double *ly)
{
    if((*z)++ % 2 == 0 && numpop > 1)
    {
        *lx = 350;
        *ly = page_height - 190;
    }
    else
    {
        *lx = 100;
        page_height -= 180;
        if(page_height-50 < 180)
        {
            pdf_new_page("");
            page_height = (double) HPDF_Page_GetHeight(page) - 10.0;
        }
        *ly = page_height - 190.0;
    }
}

///
/// plot posterior distribution, returns page_height
/// 
double pdf_loci_histogram(world_fmt *world)
{
    long locus         = world->loci > 1 ? world->loci : 0;
    return pdf_locus_histogram(world,locus);
}

/// plots posterior for a single locus
double pdf_locus_histogram(world_fmt *world, long locus)
{
    bayes_fmt * bayes  = world->bayes;
    const long numpop  = world->numpop;
    const long numpop2 = world->numpop2;
    const long np2x    = numpop2 + (world->bayes->mu);
    const long np2xx   = np2x + 2 * world->species_model_size;
    const long np2xxx  = np2xx + world->grownum;
    const long numparam= np2xxx;
    long z             = 0;
    long numbins       = 0;
    long numbinsall    = 0;
    //long p99bins       = 0;
    long *bins         = bayes->histogram[locus].bins;
    long frompop;
    long topop;
    long pa0, pa=0;
    long rpa;
    //long i;
    double themin       = 0;
    double themax       = 0;
    double thetamin     = (double) HUGE;
    double migmin       = (double) HUGE;
    double ratemin      = (double) HUGE;
    double specmin      = (double) HUGE;
    double thetamax     = 0;
    double migmax       = 0;
    double ratemax      = 0;
    double specmax      = 0;
    double growthmin    = (double) HUGE;
    double growthmax    = (double) -HUGE; 
    double mlfmin    = (double) 0.4;
    double mlfmax    = (double) 1.0; 
    page_height  = (double) HPDF_Page_GetHeight(page);
    double page_width   = (double) HPDF_Page_GetWidth(page);
    double lx;
    double ly;
    
    char *set50        = bayes->histogram[locus].set50;
    char *set95        = bayes->histogram[locus].set95;
    char title[100];
    MYREAL meanmu      = 0.0;
    MYREAL mu          = 0.0;
    long lmu         = 0;
    long l;
    double delta;
    double *results     = bayes->histogram[locus].results;
    MYREAL *mini       = bayes->histogram[locus].minima;
    MYREAL *maxi       = bayes->histogram[locus].maxima;
    
    // set the title of the section
    if (locus < world->loci)
      sprintf(title,"%s %li","Bayesian Analysis: Posterior distribution for locus",locus + 1);
    else
      sprintf(title,"%s","Bayesian Analysis: Posterior distribution over all loci");

    pdf_title(title, page_width);
    
    switch(bayes->prettyhist)
    {
    case PRETTY_MAX:
    case PRETTY_P99:
    case PRETTY_P100:
      break;
    default:
      find_posterior_min_max(&thetamin, &thetamax, 0, numpop, bayes, locus);
      find_posterior_min_max(&migmin, &migmax, numpop, numpop2, bayes, locus);
      if (world->has_speciation)
	find_posterior_min_max(&specmin, &specmax, np2x, np2xx, bayes, locus);
      if(world->has_growth)
	find_posterior_min_max(&growthmin, &growthmax, np2xx, np2xxx, bayes, locus);
      if(bayes->mu)
	{
	  find_posterior_min_max(&ratemin, &ratemax, numpop2, np2x, bayes, locus);
	}
    }
    if (bayes->mu && locus == world->loci)
      { 
	if(world->loci>1)
	  {
	    meanmu = 0.;
	    for(l=0;l<world->loci;l++)
	      {
		meanmu += world->options->meanmu[l];
	      }
	    meanmu /= world->loci;
	    if (meanmu < SMALL_VALUE)
	      meanmu = 10e-8;
	    lmu = (long) (floor( log10(bayes->histogram[locus].modes[numpop2] * meanmu)));
	    mu = meanmu * pow(10. , -lmu);
	  }
	else
	  {
	    lmu = (long) (floor( log10(bayes->histogram[locus].modes[numpop2] * world->options->meanmu[0])));
	    mu = world->options->meanmu[0] * pow(10., -lmu);
	  }
	ratemin = (double) (ratemin * mu);
	ratemax = (double) (ratemax * mu);
      }
    lx = 100;
    ly = page_height - 190;
    numbinsall = 0;
    for(pa0=0; pa0 < numparam; pa0++)
    {
      if(shortcut(pa0,world,&pa))
	{
	  continue;
	}
      if(pa < pa0)
	continue; // does not print multiple copies for 'm' and 's'
      rpa=pa;
      numbins = numbinsall;
      numbinsall += bins[rpa];
      if(rpa<numpop)
	    {
	      pdf_print_contents_at(lx-30,ly+125, "Freq");
	      symbol_Theta(lx+80, ly-25, 12, pa0+1);//frompop+1);
	      themin = thetamin;
	      themax = thetamax;
	    }
      else if (rpa < numpop2)
	    {
	      m2mm(pa0, numpop, &frompop, &topop);
	      pdf_print_contents_at(lx-30,ly+125, "Freq");
	      symbol_M(lx+80, ly-25, 12, frompop+1, topop+1, world->options->usem);
	      themin = migmin;
	      themax = migmax;
            }
      else if (bayes->mu && rpa==numpop2)
            {
	      pdf_print_contents_at(lx-30,ly+125, "Freq");
	      symbol_R(lx+80, ly-25, 12, lmu);
	      themin = ratemin;
	      themax = ratemax;
            }
      else if (world->has_speciation && rpa < np2xx)
	    {
	      species_fmt * s = get_which_species_model(pa,world->species_model,world->species_model_size);
	      frompop = s->from;
	      topop   = s->to;
	      if(pa == s->paramindex_mu)
		{
		  themin = (double) s->min;
		  themax = (double) s->max;
		  pdf_print_contents_at(lx-30,ly+125, "Freq");
		  symbol_D(lx+80, ly-25, 12, frompop+1, topop+1);
		}
	      else
		{
		  pdf_print_contents_at(lx-30,ly+125, "Freq");
		  symbol_S(lx+80, ly-25, 12, frompop+1, topop+1);
		  themin = (double) s->sigmamin;
		  themax = (double) s->sigmamax;
		}	
	    }
      else if (world->has_growth && rpa < np2xxx)
	{
	  long d=rpa-np2xx;
	  while (d < world->options->growpops_numalloc && world->options->growpops[d]==0)
	    d++;
	  pdf_print_contents_at(lx-30,ly+125, "Freq");
	  symbol_Growth(lx+80, ly-25, 12, d+1);
	  themin = growthmin;
	  themax = growthmax;
	}
    	delta = (double) (maxi[rpa] - mini[rpa])/bins[rpa];
	if (bayes->mu && rpa == numpop2)
	  delta = (double) (delta * mu);
	pdf_pretty_histogram(pa, rpa, numbins, results, 0, set50,
			     set95, bins, delta, mini, maxi,
			     lx, ly, world, themin, themax);
	pdf_advance_hist(world->numpop, &z, &lx, &ly);
    }
    return page_height;
}

///
/// prints right-adjusted text, based on margin-width x and stay on the line
void pdf_printf_right(double x, double y, char string[],...)
{
    //    double     page_width = HPDF_Page_GetWidth(page);
    double w;
    char message[LINESIZE];
	char fp[LINESIZE];
	va_list args;
    va_start (args, string);
    vsprintf (message, string, args);
    va_end (args);
    sprintf(fp,"%s",message);
    w = (double) HPDF_Page_TextWidth(page, fp);
    pdf_print_contents_at(/*page_width-*/x-w, y, fp);
}


///
/// prints right-adjusted text, based on margin-width x and jumps to next line
void pdf_printf_right_next(double x, double *y, char string[],...)
{
  double page_width = (double) HPDF_Page_GetWidth(page);
    double w;
    char message[LINESIZE];
    char fp[LINESIZE];
    va_list args;
    va_start (args, string);
    vsprintf (message, string, args);
    va_end (args);
    sprintf(fp,"%s",message);
    w = (double) HPDF_Page_TextWidth(page, fp);
    pdf_print_contents_at(page_width-x-w, *y, fp);
    pdf_advance(y);
}

///
/// prints right-aligned mutable text. Like vprintf
void pdf_printf_ralign(double rx, double y, char string[], ...)
{
    double w;
    char message[LINESIZE];
	char fp[LINESIZE];
	va_list args;
    va_start (args, string);
    vsprintf (message, string, args);
    va_end (args);
	sprintf(fp,"%s",message);
    w = (double) HPDF_Page_TextWidth(page, fp);
    pdf_print_contents_at(rx-w, y, fp);
}



///
/// Prints mutable char; like putc.
void pdf_putc(double *x, double *y, double leftborder, double rightborder, char message)
{
    char fp[100];
    double w;
    
    sprintf(fp,"%c",message);
    fp[1]='\0';
    w = (double) HPDF_Page_TextWidth(page, fp);
    if(*x + w >= rightborder)
    {
        pdf_advance(y);
        *x = leftborder;
    }
    pdf_print_contents_at(*x, *y, fp);
    *x += w;
}

///
/// Prints aligned text mutable text; like vprintf.
/// The align parameter defines the alignment of the coordinate
/// when align = 'L' it is left aligned, 'C' center aligned, and 'R' rightaligned
void pdf_printf(double x, double y, char align, char string[], ...)
{
    
    char message[LINESIZE];
    char fp[LINESIZE];
    double w;
    va_list args;
    
    va_start (args, string);
    vsprintf (message, string, args);
    va_end (args);
    
    sprintf(fp,"%s",message);
    
    w = (double) HPDF_Page_TextWidth(page, fp);
    switch(align)
    {
        case 'C':
            w /= 2.;
            break;
        case 'R':
            break;
        case 'L':
        default:
            w = 0.;
            break;
    }
    pdf_print_contents_at(x+w, y, fp);
}


///
/// prints mutable text and advances to next line.
void pdf_printf_next(double x, double *y, char string[], ...)
{
    char message[LINESIZE];
	char fp[LINESIZE];
	va_list args;
    va_start (args, string);
    vsprintf (message, string, args);
    va_end (args);
	sprintf(fp,"%s",message);
    pdf_print_contents_at(x, *y, fp);
    pdf_advance(y);
}

///
/// prints mutable text at x,y and changes x and y when reaching end of line
void pdf_printf_cell(double *x, double *y, double width, char string[], ...)
{
    char message[LINESIZE];
    char fp[LINESIZE];
    double w;
    //double ww;
    va_list args;
    
    va_start (args, string);
    vsprintf (message, string, args);
    va_end (args);
    sprintf(fp,"%s",message);
    w = (double) HPDF_Page_TextWidth(page, fp);
    if((w + *x) >=  width)
    {
        *x = 55;
        pdf_advance(y);
    }
    pdf_print_contents_at(*x, *y, fp);
    *x += w+5;
}

void    pdf_print_connection_table (world_fmt * world, option_fmt *options, data_fmt * data)
{
    const double offset = 20;
    long i;
    long j;
    long z;
    double right_margin = (double) HPDF_Page_GetWidth(page) - 2.0*55.0 - offset;
    left_margin = 55;
    double left_margin2 = 60;
    double lx = left_margin + 105;
    pdf_advance(&page_height);
    pdf_printf_next(left_margin, &page_height,   "Connection matrix:");
    pdf_printf_next (left_margin2, &page_height, "m = average (average over a group of Thetas or M,");
    pdf_printf_next (left_margin2, &page_height, "s = symmetric migration M, S = symmetric 4Nm,");
    pdf_printf_next (left_margin2, &page_height, "0 = zero, and not estimated,");
    pdf_printf_next (left_margin2, &page_height, "* = migration free to vary, Thetas are on diagonal");
    pdf_printf_next (left_margin2, &page_height, "d = row population split off column population, D = split and then migration");

    pdf_advance(&page_height);
    pdf_printf (left_margin, page_height, 'L', "Population");
    for (i = 0; i < data->numpop; i++)
    {
        pdf_printf (lx, page_height, 'L', "%3li", options->newpops[i]);
        lx += offset;
        if(lx > right_margin)
        {
            lx = left_margin + 110;
            pdf_advance(&page_height);
        }
    }
    //xcode lx = left_margin + 110;
    for (i = 0; i < data->numpop; i++)
    {
        lx = left_margin + 110;
        pdf_advance(&page_height);
        pdf_printf (left_margin, page_height, 'L', "%3li %-15.15s",options->newpops[i],  data->popnames[i]);
        for (j = 0; j < data->numpop; j++)
        {
	  //z = mm2m(options->newpops[j]-1,options->newpops[i]-1,world->numpop);
	  z = world->numpop * (options->newpops[i]-1) + (options->newpops[j]-1);	  
	  pdf_printf(lx, page_height, 'L', " %c", world->options->custm[z]);
	  lx += offset;
	  if(lx > right_margin)
            {
	      lx = left_margin + 90;
	      pdf_advance(&page_height);
            }
        }        
    }
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    if (world->has_growth) //!(world->options->growpops[0] == 0 && world->options->growpops_numalloc==1))
      {
	for (i = 0; i < data->numpop; i++)
	  {
	    lx = left_margin + 110;
	    pdf_advance(&page_height);
	    pdf_printf (left_margin, page_height, 'L', "%3li %-15.15s",options->newpops[i],  data->popnames[i]);
	    pdf_printf(lx, page_height, 'L', "%s", (world->options->growpops[options->newpops[i]-1]>0) ?  "growing/shrinking" : "constant size");
	  }
      }
    pdf_advance(&page_height);
    pdf_advance(&page_height);
}


///
/// prints order of parameters for PDF output file
void pdf_print_param_order(world_fmt *world)
{
    left_margin = 60;
    long pa;
    long pa0;
    long z=1;
    long numpop = world->numpop;
    long numpop2 = world->numpop2;
    long frompop, topop;
    long frompop2, topop2;
    char *custm2 = world->options->custm2;
    boolean usem = world->options->usem;
    //bayes_fmt *bayes = world->bayes;
    long numparam = world->numpop2 + world->bayes->mu + world->species_model_size * 2;
    pdf_advance(&page_height);
    pdf_printf_next(left_margin, &page_height,"Order of parameters:");
    //  fprintf(bayesfile,"# Parameter-number Parameter\n");
    for(pa0=0;pa0<numparam+world->grownum;pa0++)
    {
      if(shortcut(pa0,world,&pa))
	{
	  continue;
	}
      if(pa0 < numpop)
	{
	  if((pa0 == pa) && (custm2[pa0] == '*'))
	    {
	      pdf_printf(left_margin, page_height,'L',"%4li ",z++);
	      symbol_Theta(left_margin+80, page_height,12,pa0+1);
	      pdf_printf(left_margin+230, page_height,'L',"<displayed>");
	    }
	  else  /*theta option other than **/
	    {
	      pdf_printf(left_margin, page_height,'L',"%4li ",z++);
	      symbol_Theta(left_margin+80, page_height,12,pa0+1);
	      pdf_printf(left_margin+120, page_height,'L'," = ");
	      symbol_Theta(left_margin+150, page_height,12,pa+1);
	      pdf_printf(left_margin+190, page_height,'L',"[%c]",custm2[pa0]);
	      if(pa==pa0)
		pdf_printf(left_margin+230, page_height,'L',"<displayed>");
	    }
	}
      else /*migration or rate option or speciation or growth*/
	{
	  // do we estimate mutation rate changes?
	  if(world->bayes->mu && pa0==numpop2)
	    {
	      pdf_printf(left_margin, page_height,'L',"%4li ",z++);
	      pdf_printf(left_margin+80, page_height,'L', "Rate");
	      pdf_printf(left_margin+230, page_height,'L',"<displayed>");
	    }
	  else /*migration*/
	    {
	      if(pa0<numpop2)
		{
		  m2mm(pa0,numpop,&frompop,&topop);
		  if((pa0==pa) && (custm2[pa0]=='*'))
		    {
		      if(usem)
			{
			  pdf_printf(left_margin, page_height,'L',"%4li ",z++);
			  symbol_M(left_margin+80, page_height,12,frompop+1,topop+1,TRUE);
			  pdf_printf(left_margin+230, page_height,'L',"<displayed>");
			}
		      else
			{
			  pdf_printf(left_margin, page_height,'L',"%4li ",z++);
			  symbol_M(left_margin+80, page_height,12,frompop+1,topop+1,FALSE);
			  pdf_printf(left_margin+120, page_height,'L'," = ");
			  symbol_Theta(left_margin+150, page_height,12,topop+1);
			  symbol_M(left_margin+190, page_height,12,frompop+1,topop+1,TRUE);
			  pdf_printf(left_margin+230, page_height,'L',"<displayed>");
			}
		    }
		  else /* migration option other than * */
		    {
		      m2mm(pa,numpop,&frompop2,&topop2);
		      if(usem)
			{
			  pdf_printf(left_margin, page_height,'L',"%4li ",z++);
			  symbol_M(left_margin+80, page_height,12,frompop+1,topop+1,TRUE);
			  pdf_printf(left_margin+120, page_height,'L'," = ");
			  symbol_M(left_margin+150, page_height,12,frompop2+1,topop2+1,TRUE);
			  pdf_printf(left_margin+190, page_height,'L',"[%c]",custm2[pa0]);
			  if(pa==pa0)
			    pdf_printf(left_margin+230, page_height,'L',"<displayed>");
			}
		      else
			{
			  symbol_M(left_margin+80, page_height,12,frompop+1,topop+1,FALSE);
			  pdf_printf(left_margin+120, page_height,'L'," = ");
			  symbol_Theta(left_margin+150, page_height,12,topop+1);
			  symbol_M(left_margin+190, page_height,12,frompop+1,topop+1,TRUE);
			  pdf_printf(left_margin+240, page_height,'L'," = ");
			  symbol_M(left_margin+290, page_height,12,frompop2+1,topop2+1,FALSE);
			  pdf_printf(left_margin+330, page_height,'L'," = ");
			  symbol_Theta(left_margin+350, page_height,12,topop2+1);
			  symbol_M(left_margin+390, page_height,12,frompop2+1,topop2+1,TRUE);
			  pdf_printf(left_margin+430, page_height,'L',"[%c]",custm2[pa0]);
			  if(pa==pa0)
			    pdf_printf(left_margin+460, page_height,'L',"<displayed>");
			  
			}
		    }
		}
	      else  /*has speciation or growth*/
		{
		  if (world->has_speciation && pa < numparam)
		    {
		      species_fmt * s = get_which_species_model(pa, world->species_model, world->species_model_size);
		      frompop = s->from;
		      topop = s->to;
		      pdf_printf(left_margin, page_height,'L',"%4li ",z++);
		      if(pa == s->paramindex_mu)
			{
			  symbol_D(left_margin+80, page_height,12,frompop+1,topop+1);
			  pdf_printf(left_margin+230, page_height,'L',"<displayed>");
			}
		      else
			{
			  symbol_S(left_margin+80, page_height,12,frompop+1,topop+1);
			  pdf_printf(left_margin+230, page_height,'L',"<displayed>");
			}
		      pdf_advance(&page_height);
		    }
		  else
		    {
		      if (world->has_growth && pa >= numparam)
			{
			  pdf_printf(left_margin, page_height,'L',"%4li ",z++);
			  long d=pa-numparam;
			  while (d < world->options->growpops_numalloc && world->options->growpops[d]==0)
			    d++;
			  symbol_Growth(left_margin+80, page_height, 12, d+1);
			  pdf_printf(left_margin+230, page_height,'L',"<displayed>");
			}

		      //?????? error("problem with PDF histogram plotting");
		    }
		}
	    }
	}
      pdf_advance(&page_height);
    }
    pdf_advance(&page_height);
    pdf_advance(&page_height);
}



void    pdf_print_distance_table (world_fmt * world, option_fmt * options, data_fmt * data)
{
  (void) world;
    const double offset = 60;
    long i;
    long j;
    double right_margin = (double) HPDF_Page_GetWidth(page) - 2.0*55.0 - offset;
    left_margin = 55;
    double lx=left_margin + 80;
    
    if(!options->geo)
        return;
    
    pdf_advance(&page_height);
    pdf_printf_next(left_margin, &page_height, "Distance among populations:");
    pdf_advance(&page_height);
    pdf_printf (left_margin, page_height, 'L', "Population");
    for (i = 0; i < data->numpop; i++)
    {
        pdf_printf (lx, page_height, 'L', "%3li", options->newpops[i]);
        lx += offset;
        if(lx > right_margin)
        {
            lx = left_margin + 80;
            pdf_advance(&page_height);
        }
    }
    for (i = 0; i < data->numpop; i++)
    {
        lx = left_margin + 80;
        pdf_advance(&page_height);
        pdf_printf (left_margin, page_height, 'L', "%3li %-15.15s",options->newpops[i],  data->popnames[i]);
        for (j = 0; j < data->numpop; j++)
        {
            pdf_printf(lx, page_height, 'L', " %10.4f ", data->ogeo[j][i]);
            lx += offset;
            if(lx > right_margin)
            {
                lx = left_margin + 80;
                pdf_advance(&page_height);
            }
        }
    }
    pdf_advance(&page_height);
    pdf_advance(&page_height);
}



void pdf_print_options(world_fmt * world, option_fmt *options, data_fmt * data)
{
  //const char text[8][20]={"All", "Multiplier", "Exponent.", "Exp window", "Gamma", "Uniform", "Truncated Normal", "-"};
    char ptypename[100];
    double page_width;
    //page_height = *orig_page_height;
    left_margin = 55;
    double right_margin = 300.;
    double w;
    //double width[11]={0, 85, 100, 105, 150, 210, 260, 300, 360, 420, 480};
    double width[11]={0, -20, 100, 110, 180, 230, 290, 340, 400, 440, 500};
    double width2[2]={0, 240};
    char *title = "Options";
    long i, j, tt, ii;
    //long pop;
    char mytext[LINESIZE];
    char mytext1[LINESIZE];
    char mytext2[LINESIZE];
    char mytext3[LINESIZE];
    char mytext4[LINESIZE];
    char mytext5[LINESIZE];
    char mytext6[LINESIZE];
    char seedgen[LINESIZE], spacer[LINESIZE];
    //char fromstring[LINESIZE];
    //char tostring[LINESIZE];
    //char **paramtgen;
    boolean inheritance_table = FALSE;
    char * priorkind = (char *) mycalloc(LINESIZE, sizeof(char));
    //xcode page_width = HPDF_Page_GetWidth(page) - 120;
    //charvec2d(&paramtgen,PRIOR_SIZE,LINESIZE);
    char *paramtgen, *parammgen;
    paramtgen = (char *) mycalloc(2*LINESIZE,sizeof(char));
    parammgen = paramtgen + LINESIZE;

    if (options->datatype != 'g')
    {
        switch ((short) options->autoseed)
        {
            case AUTO:
                strcpy (seedgen, "with internal timer");
                strcpy (spacer, "  ");
                break;
            case NOAUTOSELF:
                strcpy (seedgen, "from parmfile");
                strcpy (spacer, "      ");
                break;
            case NOAUTO:
                strcpy (seedgen, "from seedfile");
                strcpy (spacer, "      ");
                break;
            default:
                strcpy (seedgen, "ERROR");
                strcpy (spacer, " ");
                break;
        }
	fill_printvar_startparam(options, &paramtgen, &parammgen);
    }
    HPDF_Page_SetFontAndSize (page, helvob, 18.0);
    w = (double) HPDF_Page_TextWidth(page, title);
    page_width = (double) HPDF_Page_GetWidth(page);
    right_margin = page_width - 55.0;
    HPDF_Page_SetRGBStroke(page, 0.f, 0.f, 0.f);
    HPDF_Page_SetRGBFill(page, 0.f, 0.f, 0.f);
    HPDF_Page_SetLineWidth(page, 1);
    pdf_print_contents_at((page_width - w)/2, page_height, title);
    page_height -= 20;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    page_height -= 20;
    HPDF_Page_SetFontAndSize (page, helv, 10.);
    
    switch (options->datatype)
    {
        case 'a':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            pdf_printf_right_next(left_margin, &page_height,"Allelic data");
            pdf_print_contents_at(left_margin, page_height,"Missing data:");
            pdf_printf_right_next(left_margin, &page_height,"%s\n",options->include_unknown ? "included" : "not included");
            break;
        case 'b':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            pdf_printf_right_next(left_margin, &page_height,"Microsatellite data [Brownian motion]\n");
            pdf_print_contents_at(left_margin, page_height,"Missing data:");
            pdf_printf_right_next(left_margin, &page_height,"%s\n",options->include_unknown ? "included" : "not included");
            break;
        case 'm':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            if(options->msat_option == SINGLESTEP)
                pdf_printf_right_next(left_margin, &page_height,"Microsatellite data [Singlestep model]\n");
            else
                pdf_printf_right_next(left_margin, &page_height,"Microsatellite data [Multistep model (Tune=%f, P_increase=%f)]\n",
                                      options->msat_tuning[0], options->msat_tuning[1]);
            pdf_print_contents_at(left_margin, page_height,"Missing data:");
            pdf_printf_right_next(left_margin, &page_height,"%s\n",options->include_unknown ? "included" : "not included");
            break;
        case 's':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            pdf_printf_right_next(left_margin, &page_height,"DNA sequence data\n");
            break;
        case 'n':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            pdf_printf_right_next(left_margin, &page_height,"Single nucleotide polymorphism data\n");
            break;
        case 'h':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            pdf_printf_right_next(left_margin, &page_height,"Single nucleotide polymorphism data(Hapmap formatting)\n");
            break;
        case 'u':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            pdf_printf_right_next(left_margin, &page_height,"Single nucleotide polymorphism data (PANEL)\n");
            break;
        case 'f':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            pdf_printf_right_next(left_margin, &page_height,"Ancestral state method\n");
            break;
        case 'g':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            pdf_printf_right_next(left_margin, &page_height,"Genealogy summary of an older run\n");
            break;
    }
    
    pdf_advance(&page_height);
    pdf_print_contents_at(left_margin, page_height,"Inheritance multipliers in use for Thetas:");
    pdf_advance(&page_height);
    //left_margin += HPDF_Page_TextWidth(page, "Inheritance scalers in use for Thetas: ");
    for(i=0;i<data->loci;i++)
    {
        if(options->inheritance_scalars[i]<1.0 || options->inheritance_scalars[i]>1.0)
        {
            inheritance_table = TRUE;
            break;
        }
    }
    if (inheritance_table)
    {
        for(i=0;i<data->loci;i++)
        {
            if(i<options->inheritance_scalars_numalloc)
	      pdf_printf_cell(&left_margin, &page_height, page_width-55.0, "%2.2f", (double) options->inheritance_scalars[i]);
            else
	      pdf_printf_cell(&left_margin, &page_height, page_width-55.0, "%2.2f", (double) options->inheritance_scalars[options->inheritance_scalars_numalloc-1]);
        }
    }
    else
    {
        pdf_print_contents_at(left_margin, page_height,"All loci use an inheritance multiplier of 1.0");
    }
    pdf_advance(&page_height);
    //left_margin = *orig_left_margin;
    //pdf_print_contents_at(left_margin, page_height,"[The locus with a scaler of 1.0 used as reference]");
    //pdf_advance(&page_height);
    
    
    if(options->randomsubset > 0)
    {
        pdf_advance(&page_height);
        pdf_print_contents_at(left_margin, page_height,"Data set was subsampled: used a random sample of size:");
        pdf_advance(&page_height);
        if (options->randomsubsetseed > 0)
            pdf_printf_right_next(left_margin, &page_height,"%5li and seed %8li", options->randomsubset,options->randomsubsetseed);
        else
            pdf_printf_right_next(left_margin, &page_height,"%5li", options->randomsubset);
        pdf_advance(&page_height);
    }
    
    if (options->datatype != 'g')
    {
        pdf_print_contents_at(left_margin, page_height,"Random number seed:");
        pdf_printf_right_next(left_margin, &page_height,"(%s)%s%20li", seedgen, " ",
                              options->saveseed);
        pdf_printf_next(left_margin, &page_height,"Start parameters:");
        pdf_advance(&page_height);
        pdf_print_contents_at(left_margin, page_height,"Theta values were generated");
        pdf_printf_right_next(left_margin, &page_height," %s", paramtgen);
        if (options->startguess[THETAPRIOR][0] == OWN)
        {
            //xcode ii = 0;
            left_margin += 10;
            pdf_print_contents_at(left_margin, page_height,"Theta = ");
            left_margin += (double) HPDF_Page_TextWidth(page, "Theta = ");
            for (i = 0; i < options->startparam.numtheta; i++)
            {
	      pdf_printf_cell(&left_margin, &page_height, page_width-55.0, "%.5f", (double) options->startparam.theta[i]);
            }
        }
        
        //left_margin = *orig_left_margin;
        pdf_advance(&page_height);
        
        if(options->usem)
            pdf_print_contents_at(left_margin, page_height,"M values were generated");
        else
            pdf_print_contents_at(left_margin, page_height,"xNm values were generated");
        
        pdf_printf_right_next(left_margin, &page_height," %s", parammgen);
        
        //left_margin = *orig_left_margin + 10;
        if (options->startguess[MIGPRIOR][0] == OWN)
        {
            tt = 0;
            if (options->usem)
                pdf_print_contents_at(left_margin, page_height,"M-matrix:");
            else
                pdf_print_contents_at(left_margin, page_height,"xNm-matrix:");
            pdf_advance(&page_height);
            if (options->startparam.nummig == 1)
            {
                pdf_printf_cell(&left_margin, &page_height, page_width-55.0,
                                "%5.2f [all are the same]", (double) options->startparam.mig[tt]);//xcode replaced tt++ with tt
            }
            else
            {
                for (i = 0; i < world->numpop; i++)
                {
                    //xcode  ii = 0;
                    for (j = 0; j < world->numpop; j++)
                    {
                        if (i != j)
                        {
                            pdf_printf_cell(&left_margin, &page_height, page_width-55.0,
                                            " %5.*f, ",options->usem ? 1 : 5, (double) options->startparam.mig[tt++]);
                        }
                        else
                        {
                            pdf_printf_cell(&left_margin,
                                            &page_height, page_width-55.0, "   -   ");
                        }
                    }
                    //left_margin = *orig_left_margin;
                    pdf_advance(&page_height);
                }
                //left_margin = *orig_left_margin;
                pdf_advance(&page_height);
            }
        }
    }
    pdf_print_connection_table (world, options, data);
    pdf_print_distance_table (world, options, data);
    pdf_print_param_order(world);
    //left_margin = *orig_left_margin;
    if (options->gamma)
    {
        pdf_print_contents_at(left_margin, page_height,"Mutation rate among loci:");
        pdf_printf_right_next(left_margin, &page_height,"from a Gamma distribution\n");
        pdf_printf_right_next(left_margin, &page_height,"Initial scale parameter alpha = %f\n",
                              options->alphavalue);
        if (options->custm[world->numpop2] == 'c')
        {
            pdf_printf_next(left_margin + 300, &page_height,"and is constant [will not be estimated]\n");
        }
    }
    else
    {
        if(options->bayesmurates)
        {
            pdf_printf_right_next(left_margin, &page_height,"Mutation rate is estimated %s", world->loci > 1 ?
                                  "for all loci" : "");
        }
        else
        {
            pdf_print_contents_at(left_margin, page_height,"Mutation rate among loci:");
            if (options->murates && world->loci > 1)
            {
                if(options->murates_fromdata)
                    pdf_printf_right_next(left_margin, &page_height,"Varying ([crudely] estimated from data)");
                else
                    pdf_printf_right_next(left_margin, &page_height,"Varying (user input)");
                pdf_print_contents_at(left_margin + 10, page_height,"Rates per locus: ");
                ii=0;
                for (i = 0; i < world->loci-1; i++)
                {
                    sprintf(mytext,"%.5f, ", options->mu_rates[i]);
                    if (i % 6 == 5)
                    {
                        ii=0;
                        pdf_advance(&page_height);
                        pdf_print_contents_at(left_margin + 100 + ii * 60, page_height,mytext);
                    }
                    else
                    {
                        pdf_print_contents_at(left_margin + 100 + ii * 60, page_height,mytext);
                    }
                    ii++;
                }
                if(i % 6 == 5)
                {
                    ii=0;
                    pdf_advance(&page_height);
                }
                pdf_printf_next(left_margin + 100 + ii * 60, &page_height,"%.5f", options->mu_rates[i]);
            }
            else
                pdf_printf_right_next(left_margin, &page_height,"Mutation rate is constant %s", world->loci > 1 ?
                                      "for all loci" : "");
        }
    }
#ifdef UEP
    if (options->uep)
    {
        pdf_printf_next(left_margin, &page_height,"0/1 polymorphism analysis, with 0/1 data in file:");
        pdf_printf_right_next(left_margin, &page_height,"%s\n", options->uepfilename);
        pdf_printf_right_next(left_margin, &page_height,"with forward mutation rate %f*mu\n",
                              options->uepmu);
        pdf_printf_right_next(left_margin, &page_height,"with back mutation rate %f*mu\n",
                              options->uepnu);
        pdf_printf_right_next(left_margin, &page_height,"with base frequencies \"0\"=%f and \"1\"=%f\n",
                              options->uepfreq0,options->uepfreq1);
    }
#endif
    pdf_advance(&page_height);
    
    pdf_print_contents_at(left_margin, page_height,"Analysis strategy:");
    pdf_printf_right_next(left_margin, &page_height,"Bayesian inference");
    //pdf_advance(&page_height);
    pdf_print_contents_at(left_margin, page_height," -Population size estimation:");
    if (options->mlalpha<1.0)
      sprintf(mytext6,"Mittag-Leffler with alpha=%.2f",options->mlalpha);
    else
      sprintf(mytext6,"Exponential Distribution");
    pdf_printf_right_next(left_margin, &page_height,"%s",mytext6);
    //pdf_advance(&page_height);

    if(world->has_growth)
      {
	pdf_print_contents_at(left_margin, page_height," -Population growth estimation: Growth [Exponential]\n");
	pdf_advance(&page_height);
      }
        if (world->numpop>1)
      {
	if(world->has_migration)
	  {
	    pdf_print_contents_at(left_margin, page_height," -Geneflow estimation:");
	    pdf_printf_right_next(left_margin, &page_height,"%s", mytext6);
	    //pdf_advance(&page_height);
	  }
	if(world->has_speciation)
	  {
	    switch(world->species_model_dist)
	      {
	      case WEIBULL_DIST:
		pdf_print_contents_at(left_margin, page_height," -Divergence time estimation:");
		pdf_printf_right_next(left_margin, &page_height,"Weibull Distribution (mean and spread parameter)");
		//pdf_advance(&page_height);
		break;
	      case EXP_DIST:
		pdf_print_contents_at(left_margin, page_height," -Divergence time estimation:");
		pdf_printf_right_next(left_margin, &page_height,"Exponential Distribution (mean)");
		break;
	      case NORMAL_DIST:
	      case NORMALSHORTCUT_DIST:
	      default:
		pdf_print_contents_at(left_margin, page_height," -Divergence time estimation:");
		pdf_printf_right_next(left_margin, &page_height,"Normal Distribution (mean and standard dev.)");
		//pdf_advance(&page_height);
		break;
	      }
	  }
      }	
    pdf_advance(&page_height);
    pdf_print_contents_at(left_margin, page_height,"Proposal distributions for parameter");
    pdf_advance(&page_height);
    pdf_print_tableline(width2, "%s %s", "Parameter", "Proposal");
    pdf_advance(&page_height);
    pdf_print_tableline(width2, "%s %s", "Theta",
			is_proposaltype(options->slice_sampling[THETAPRIOR]));
    pdf_advance(&page_height);
    if (world->has_migration)
      {
	pdf_advance(&page_height);
	pdf_print_tableline(width2, "%s %s", options->usem ? "M" : "xNm",
			    is_proposaltype(options->slice_sampling[MIGPRIOR]));
      }
    if (world->has_speciation)
      {
	pdf_advance(&page_height);
	pdf_print_tableline(width2, "%s %s", "Divergence",
			    is_proposaltype(options->slice_sampling[SPLITPRIOR]));
	pdf_advance(&page_height);
	pdf_print_tableline(width2, "%s %s", "Divergence Spread",
			    is_proposaltype(options->slice_sampling[SPLITSTDPRIOR]));
      }
    if (world->has_growth)
      {
	pdf_advance(&page_height);
	pdf_print_tableline(width2, "%s %s", "Growth",
			    is_proposaltype(options->slice_sampling[GROWTHPRIOR]));
      }
    if(options->bayesmurates)
      {
	pdf_advance(&page_height);
	pdf_print_tableline(width2, "%s %s", "Rate",
			    is_proposaltype(options->slice_sampling[RATEPRIOR]));
      }
    pdf_advance(&page_height);
    pdf_print_tableline(width2, "%s %s", "Genealogy",
			"Metropolis-Hastings");
    pdf_advance(&page_height);
    pdf_advance(&page_height);    
    pdf_print_contents_at(left_margin, page_height,"Prior distribution for parameter");
    pdf_advance(&page_height);
    pdf_print_tableline(width, "%s %s %s %s %s %s %s %s %s %s %s", "Parameter", " ", " ", " ", "Prior", "Minimum",  "Mean*",  "Maximum", "Delta", "Bins","UpdateFreq");
    pdf_advance(&page_height);

    prior_fmt *p = options->bayes_priors;
    long pnum = options->bayes_priors_num;
    long z=1;
    long numparam=0;
    //for(i=0; i < options->bayes_priors_num; i++)
    for(i=0; i < pnum; i++)
      {
	long pa=i;
    	if(shortcut(i,world,&pa))
	  continue;
	else
	  numparam++;
      }
    sprintf(mytext5,"%5.5f", (world->options->choices[1] - world->options->choices[0])/numparam);
    for(i=0; i < pnum; i++)
      {
	long pa=i;
	prior_fmt *ptr = &p[i];
	switch(ptr->type)
	  {
	  case THETAPRIOR: 	  sprintf(ptypename,"%-s","Theta");break;
	  case GROWTHPRIOR: 	  sprintf(ptypename,"%-s","Growth");break;
	  case MIGPRIOR: 	  sprintf(ptypename,"%-s",options->usem ? "M" : "xNm");break;
	  case RATEPRIOR: 	  sprintf(ptypename,"%-s","Rate modif.");break;
	  case SPECIESTIMEPRIOR: 	  sprintf(ptypename,"%-s","Splittime mean");break;
	  case SPECIESSTDPRIOR: 	  sprintf(ptypename,"%-s","Splittime std");break;
	  }
	if(shortcut(i,world,&pa))
	  continue;
	char fstr[4];
	char tstr[4];
	long from = p[pa].from;
	long to = p[pa].to;
	if (from== -1)
	  sprintf(fstr,"*  ");
	else
	  sprintf(fstr,"%2li",from+1);
	if (to == - 1)
	  sprintf(tstr,"*");
	else
	  sprintf(tstr,"%li",to+1);
	is_priorkind(ptr, priorkind);
	pdf_print_tableline(width, "%3li %s %s %s %s %s %5.5s %5.5s %9.9s %6.6s %s", z++ /*pa*/, ptypename, fstr, tstr,
			    //text[is_priortype(options->bayes_priors,options->bayes_priors_num, options->bayes_priors[pa].type)],
                           priorkind,
                           show_priormin(mytext1, ptr),
                           show_priormean(mytext2, ptr),
                           show_priormax(mytext3, ptr),
                           show_priordelta(mytext4, ptr),
                           show_priorbins(mytext, ptr),
			    mytext5
			    );
	pdf_advance(&page_height);
      }
    pdf_print_contents_at(left_margin, page_height,"[* * means priors were set globally]");   
    pdf_advance(&page_height);
    pdf_advance(&page_height);

    char smoothing1[STRSIZE] = {"--"};
    char smoothing2[STRSIZE] = {"--"};
    switch(options->smoothing[0])
      {
      case BINNING: sprintf(smoothing1,"BINNING"); break;
      case KERNEL: sprintf(smoothing1,"KERNEL SMOOTHING (window=%li)",options->smoothing_window[0]*2+1); break;
      }
    switch(options->smoothing[1])
      {
      case NONE: sprintf(smoothing2,"NO SMOOTHING"); break;
      case AVERAGE: sprintf(smoothing2,"MOVING AVERAGE SMOOTHING (window=%li)",options->smoothing_window[1]*2+1); break;
      case SAVITZKY: sprintf(smoothing2,"SAVITZKY-GOLAY SMOOTHING (window=%li)",options->smoothing_window[1]*2+1); break;
      }
    /* TODO smothing section needs better explanation and also needs to be called
       currently this is using bayes_smooth and not kernel_smooth or savitzky()*/

    pdf_print_contents_at(left_margin, page_height,"Posterior distribution:");
    pdf_advance(&page_height);
    pdf_print_contents_at(left_margin, page_height,"Parameter values were collected using MCMC, these values");
    pdf_advance(&page_height);
    sprintf(mytext,"were then used to generate the posterior histograms using %s",smoothing1);
    pdf_print_contents_at(left_margin, page_height,mytext);
    pdf_advance(&page_height);
    sprintf(mytext,"and subsequent %s for combination over loci\n", smoothing2);
    pdf_print_contents_at(left_margin, page_height, mytext);
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    //*/
    
    if (options->datatype != 'g')
    {
        pdf_print_contents_at(left_margin, page_height,"Markov chain settings:");
        pdf_printf_right_next(left_margin, &page_height, "Long chain");
        pdf_print_contents_at(left_margin, page_height,"Number of chains");
        pdf_printf_right_next(left_margin+10, &page_height,"%20li", options->lchains);
        
        pdf_print_contents_at(left_margin + 10, page_height,"Recorded steps [a]");
        pdf_printf_right_next(left_margin + 10, &page_height,"%20li", options->lsteps);
        
        pdf_print_contents_at(left_margin + 10, page_height,"Increment (record every x step [b]");
        pdf_printf_right_next(left_margin+10, &page_height,"%20li", options->lincrement);
	pdf_print_contents_at(left_margin + 10, page_height,"Number of concurrent chains (replicates) [c]");
	pdf_printf_ralign(left_margin + 340, page_height," ");
	pdf_printf_right_next(left_margin + 10, &page_height,"%20li", options->replicate ?
			      options->replicatenum : 1);
	pdf_print_contents_at(left_margin + 10, page_height,"Visited (sampled) parameter values [a*b*c]");
	// pdf_printf_ralign(left_margin + 340, page_height,"%20li",
	//                  options->ssteps * options->sincrement*(options->replicate ? options->replicatenum : 1));
	pdf_printf_right_next(left_margin + 10, &page_height,"%20li", options->lsteps * options->lincrement*(options->replicate ? options->replicatenum : 1));
        if (options->burn_in > 0)
        {
            pdf_print_contents_at(left_margin + 10, page_height,"Number of discard trees per replicate (burn-in * b)");
	    pdf_printf_right_next(left_margin + 10, &page_height,"%c%li", options->burnin_autostop, (long) options->burn_in * options->lincrement);
        }
        if (options->movingsteps)
        {
            pdf_print_contents_at(left_margin + 10, page_height,"Forcing percentage of new genealogies");
            //    pdf_printf(left_margin + 360, page_height,"%5.2f",   (MYREAL) options->acceptfreq);
            pdf_printf_right_next(left_margin + 10, &page_height,"%5.2f",   (MYREAL) options->acceptfreq);
        }
        if (options->lcepsilon < LONGCHAINEPSILON)
        {
            pdf_print_contents_at(left_margin + 10, page_height,"Forcing parameter-likelihood improvement");
            //    pdf_printf(left_margin + 360, page_height,"%20.5f",   (MYREAL) options->lcepsilon);
            pdf_printf_right_next(left_margin + 10, &page_height,"%20.5f",   (MYREAL) options->lcepsilon);
        }
        
        pdf_advance(&page_height);
        
	if(options->heating > 0)
	  pdf_printf_next(left_margin, &page_height,"Multiple Markov chains:");
               
        if (options->heating > 0)
        {
            pdf_printf(left_margin+10, page_height,'L', "%s heating scheme", options->adaptiveheat!=NOTADAPTIVE ? ( options->adaptiveheat==STANDARD ? "Adaptive_standard" : "Bounded_adaptive") : "Static");
            pdf_printf_right_next(left_margin, &page_height, "%li chains with %s temperatures",
                                  options->heated_chains, options->adaptiveheat!=NOTADAPTIVE ? "start values" : "" );
            for (i = options->heated_chains - 1; i >= 0; i--)
                pdf_printf_right(right_margin - i * 50, page_height,"%5.2f ", options->heat[i]);
            pdf_advance(&page_height);
            
            pdf_printf_right_next(left_margin, &page_height,"Swapping interval is %li\n",
                                  options->heating_interval);
        }
    }
    
    pdf_advance(&page_height);
    pdf_printf_next(left_margin, &page_height,"Print options:\n");
    if (options->datatype != 'g')
    {
        pdf_print_contents_at(left_margin + 10, page_height,"Data file:");
        pdf_printf_right_next(left_margin, &page_height,"%s", options->infilename);
	pdf_printf_right_next(left_margin, &page_height,"%s", options->parmfilename);
        pdf_print_contents_at(left_margin + 10, page_height,"Haplotyping is turned on:");
        pdf_printf_right_next(left_margin, &page_height,"%s", options->haplotyping ?
                              (options->haplotyping_report ?
                               "YES: report of haplotype probabilities"
                               : "YES: NO report of haplotype probabilities")
                              : "NO");
        pdf_print_contents_at(left_margin + 10, page_height,"Output file:");
        pdf_printf_right_next(left_margin, &page_height,"%s", options->outfilename);
        if(options->writesum)
        {
            pdf_print_contents_at(left_margin + 10, page_height,"Summary of genealogies for further run:");
            pdf_printf_right_next(left_margin, &page_height,"%s", options->sumfilename);
        }
        if(options->writelog)
        {
            pdf_print_contents_at(left_margin + 10, page_height,"Log file:");
            pdf_printf_right_next(left_margin, &page_height,"%s", options->logfilename);
        }
        
	pdf_print_contents_at(left_margin + 10, page_height,"Posterior distribution raw histogram file:");
	pdf_printf_right_next(left_margin, &page_height,"%s\n", options->bayesfilename);
	pdf_print_contents_at(left_margin + 10, page_height,"Raw data from the MCMC run:");
	pdf_printf_right_next(left_margin, &page_height,"%s\n", options->bayesmdimfilename);
        
        pdf_print_contents_at(left_margin + 10, page_height,"Print data:");
        pdf_printf_right_next(left_margin, &page_height,"%45.45s", options->printdata ? "Yes" : "No");
        
        pdf_printf(left_margin + 10, page_height,'L', "Print genealogies [only some for some data type]:");
        switch (options->treeprint)
        {
            case myNONE:
                pdf_printf_right_next(left_margin, &page_height,"None");
                break;
            case ALL:
                pdf_printf_right_next(left_margin, &page_height,"Yes, all");
                break;
            case LASTCHAIN:
                pdf_printf_right_next(left_margin, &page_height,"Yes, only those in last chain");
                break;
            case BEST:
                pdf_printf_right_next(left_margin, &page_height,"Yes, only the best");
                break;
        }
        if (options->mighist)
        {
            pdf_print_contents_at(left_margin + 10, page_height,"Histogram of the frequency of migration events");
            pdf_printf_right_next(left_margin, &page_height,"%s", options->mighistfilename);
        }
    }
    else
    {
        pdf_print_contents_at(left_margin + 10, page_height,"Data file:");
        pdf_printf_right_next(left_margin, &page_height,"%s", options->infilename);
        pdf_print_contents_at(left_margin + 10, page_height,"Output file:");
        pdf_printf_right_next(left_margin, &page_height,"%s", options->outfilename);
    }

    //*orig_page_height = page_height;
    //*orig_left_margin = left_margin;
    myfree(priorkind);
    //free_charvec2d(paramtgen);
    myfree(paramtgen);
}


void pdf_print_ratetbl (world_fmt * world, option_fmt * options, long locus, char header)
{
  (void) options;
  boolean doprint=TRUE;
  //double page_height = *orig_page_height;
  long i;
  double page_width =  (double) HPDF_Page_GetWidth(page);
  mutationmodel_fmt *s, *sp;
    long sublocus;
    const long sublocistart = world->sublocistarts[locus];
    const long sublociend   = world->sublocistarts[locus+1];
    if(header=='A' || header == 'C') 
      {
	pdf_print_contents_at(left_margin, page_height,"Site rate variation and probabilities:");
        pdf_advance(&page_height);
	pdf_printf(left_margin, page_height,'L', "Locus Sublocus Region type     Rate of change    Probability  Patch size");
        pdf_advance(&page_height);
	pdf_draw_line(50, page_height, page_width-50, page_height);
	if(header=='C')
	  {
	    pdf_advance(&page_height);
	    pdf_printf(left_margin, page_height,'L', "[compressed - only loci that are different than the one before are shown]");
	  }
      }
    for(sublocus=sublocistart; sublocus < sublociend; sublocus++)
      {
	s = &world->mutationmodels[sublocus];
	if(header=='G')
	  {
	    doprint=FALSE;
	    sp = &world->mutationmodels[sublocus-1];
	    for (i = 0; i < s->numsiterates; i++)
	      if (fabs(s->siterates[i] - sp->siterates[i]) < DBL_EPSILON)
		doprint=TRUE;	  
	  }
	if(doprint)
	  {
	    for (i = 0; i < s->numsiterates; i++)
	      {
		pdf_advance(&page_height);
		pdf_printf(left_margin, page_height,'L', "%4ld     %8ld       %9ld          %16.3f    %17.3f  %13.3f", locus + 1, sublocus + 1-sublocistart, i + 1, s->siterates[i],
			   s->siteprobs[i], 1.0 / s->lambda);
	      }
	    if (s->numsiterates > 1)
	      pdf_advance(&page_height);
	  }
      }
    //*orig_page_height = page_height;
}


void pdf_print_data_summary(world_fmt * world, option_fmt *options, data_fmt * data,
                            double *orig_page_height, double *orig_left_margin)
{
  boolean compressed=FALSE;
    long locus;
    long pop;
    long numind;
    long nummiss;
    char dstring[LINESIZE];
    char modelname[LINESIZE];
    char modelparam[LINESIZE];
    long *total;
    long *totalmiss;
    char *title = "Data summary";
    //    double w;
    //double page_height = *orig_page_height;
    //double left_margin = *orig_left_margin;
    double page_width;
    
    
    double col1 = left_margin + 300;
    double col2 = left_margin + 320;
    double col3 = left_margin + 400;
    
    
    total = (long *) mycalloc(data->loci,sizeof(long));
    totalmiss = (long *) mycalloc(data->loci,sizeof(long));
    // setup new page and title
    pdf_print_section_title(&page_width, &page_height, title);
    set_datatype_string(options->datatype, dstring);
    if(options->prioralone)
      {
	pdf_print_contents_at(left_margin, page_height,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	pdf_advance(&page_height);
	pdf_print_contents_at(left_margin, page_height,"Program is using NO DATA -- running WITHOUT DATA\n");
	pdf_advance(&page_height);
	pdf_print_contents_at(left_margin, page_height,"Option: NODATA=yes is set, to run real data remove this option\n");
	pdf_advance(&page_height);
	pdf_print_contents_at(left_margin, page_height,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	pdf_advance(&page_height);
      }
    pdf_print_contents_at(left_margin, page_height,"Data file:");
    pdf_printf_right_next(left_margin, &page_height,"%s", options->infilename);
    pdf_print_contents_at(left_margin, page_height,"Datatype:");
    pdf_printf_right_next(left_margin, &page_height,"%s", dstring);
    if(dstring[0]=='M')
      {
	if(!data->has_repeats)
	  {
	    pdf_printf_right_next(left_margin, &page_height,"[Fragment length is translated to repeats]");
	  }
	else
	  {
	    if (dstring[0]=='M')
	      pdf_printf_right_next(left_margin, &page_height, "[Data was used as repeat-length information]\n");
	  }
      }
    if(options->has_datefile)
    {
        //fprintf (file, "Sample dates:          %s\n", world->datefile);
        pdf_print_contents_at(left_margin, page_height,"Sample dates:");
        pdf_printf_right_next(left_margin, &page_height,"%s", options->datefilename);
        
        //fprintf (file, "Generations per year:  %s\n", options->generation_year);
        pdf_print_contents_at(left_margin, page_height,"Generations per year:");
        pdf_printf_right_next(left_margin, &page_height,"%f", (double) options->generation_year);
        
        //fprintf (file, "Mutationrate per year: %s", options->mutationrate_year[0]);
        pdf_print_contents_at(left_margin, page_height,"Mutationrate per year:");
        pdf_printf_right_next(left_margin, &page_height,"%f",  options->mutationrate_year[0]);
        for(locus=1; locus < options->mutationrate_year_numalloc; locus++)
        {
            pdf_printf_right_next(left_margin, &page_height,"%f",  options->mutationrate_year[locus]);
        }
    }
    pdf_print_contents_at(left_margin, page_height,"Number of loci:");
    pdf_printf_right_next(left_margin, &page_height,"%li", data->loci);
    pdf_advance(&page_height);
    pdf_print_contents_at(left_margin, page_height,"Mutationmodel:");
    pdf_advance(&page_height);
    pdf_print_contents_at(left_margin, page_height,"Locus");
    pdf_print_contents_at(left_margin+30, page_height,"Sublocus");
    pdf_print_contents_at(left_margin+90, page_height,"Mutationmodel");
    pdf_print_contents_at(left_margin+200, page_height,"Mutationmodel parameters");
    pdf_advance(&page_height);
    boolean terse = options->tersepdf;
    long printloci = data->loci;
    if(terse)
      {
	printloci = TEN;
	pdf_advance(&page_height);
	pdf_printf(left_margin, page_height,'L', "Only the first 10 loci are shown\n");
      }
    for (locus=0; locus < printloci; locus++)
      {
	long   sublocus;
	long   sublocistart = world->sublocistarts[locus];
	long   sublociend = world->sublocistarts[locus+1];
	for(sublocus=sublocistart;sublocus<sublociend;sublocus++)
	  {	  
	    get_mutationmodel_nameparam(modelname,modelparam, &world->mutationmodels[sublocus]);
	    pdf_advance(&page_height);
	    pdf_printf(left_margin, page_height,'L', "%6li     %8li           %-15.15s          %s\n",locus+1,sublocus-sublocistart+1,modelname,modelparam);
	  }
      }
    pdf_advance(&page_height);
    if(options->totalsites>0 && strchr(SNPTYPES,options->datatype))
      {
	pdf_printf(left_margin, page_height,'L', "[On average there are %li invariant sites per locus]\n",options->totalsites);
	pdf_advance(&page_height);
      }
    
    if(options->randomsubset > 0)
    {
      pdf_advance(&page_height);
      pdf_print_contents_at(left_margin, page_height,"Data set was subsampled: used a random sample of size: ");
      pdf_printf_right_next(left_margin, &page_height,"%li", options->randomsubset);
    }
    pdf_advance(&page_height);
    //print used sites per locus
    if (!(!strchr (SEQUENCETYPES, options->datatype) && options->datatype!='@'))
      {
	pdf_advance(&page_height);
	pdf_print_contents_at(left_margin, page_height,"Sites per locus");
	pdf_advance(&page_height);
	pdf_print_contents_at(left_margin, page_height,"Locus");
	pdf_print_contents_at(left_margin+100, page_height,"Sites\n");
	boolean print_siterates = FALSE;
	if (data->loci > 500)
	  compressed=TRUE;
	if (compressed)
	  {
	    long mini = 10000000;
	    long maxi = 0;
	    for(locus=0; locus< data->loci; locus++)
	      {
		long   m = 0;
		long   sublocus;
		long   sublocistart = world->sublocistarts[locus];
		long   sublociend = world->sublocistarts[locus+1];
		for(sublocus=sublocistart;sublocus<sublociend;sublocus++)
		  {
		    m = world->mutationmodels[sublocus].numsites;
		    if (world->mutationmodels[sublocus].numsiterates > 1)
		      print_siterates = TRUE;
		    if (m<mini)
		      mini = m;
		    if (m > maxi)
		      maxi = m;
		  }
	      }
	    pdf_printf(left_margin,page_height,'L',"%li loci with minimal %li and maximal %li subloci\n",
		    data->loci, mini, maxi);
	  }
	else
	  {
	    for(locus=0; locus< data->loci; locus++)
	      {
		pdf_advance(&page_height);
		pdf_printf(left_margin, page_height,'L', "%6li",locus+1);
		long z=0;
		long   sublocus;
		long   sublocistart = world->sublocistarts[locus];
		long   sublociend = world->sublocistarts[locus+1];
		for(sublocus=sublocistart;sublocus<sublociend;sublocus++)
		  {
		    if(left_margin+100+(sublocus-sublocistart)*50 > page_width)
		      {
			pdf_advance(&page_height);
			z=0;
		      }
		    if(world->mutationmodels[sublocus].numsiterates>0)
		      print_siterates = TRUE;
		    pdf_printf_ralign(left_margin+100+z*50, page_height, "%li",
				      world->mutationmodels[sublocus].numsites);
		    z++;
		  }	  
	      }
	    pdf_advance(&page_height);
	  }
	if (print_siterates)
	  {
	    pdf_advance(&page_height);
	    if(compressed)
	      pdf_print_ratetbl (world, options, 0, 'C');
	    else
	      pdf_print_ratetbl (world, options, 0, 'A');
	    for(locus=1; locus< data->loci; locus++)
	      {
		if(compressed)
		  pdf_print_ratetbl (world, options, locus, 'G');
		else
		  pdf_print_ratetbl (world, options, locus, 'F');
	      }
	    pdf_advance(&page_height);
	    pdf_advance(&page_height);
	  }
      }
    pdf_advance(&page_height);
    pdf_print_contents_at(left_margin, page_height,"Population");
    pdf_print_contents_at(col1, page_height,"Locus");
    pdf_print_contents_at(col3, page_height, "Gene copies");
    pdf_advance(&page_height);
    if (!strchr (SEQUENCETYPES, options->datatype))
    {
        pdf_print_contents_at(col3, page_height,"data");
        pdf_printf_right_next(left_margin, &page_height,"(missing)");
    }
    
    for (pop = 0; pop < data->numpop; pop++)
    {
        if (!strchr (SEQUENCETYPES, options->datatype) && options->datatype !='@')
        {
            nummiss = find_missing(data,pop,0);
            numind = data->numalleles[pop][0] - nummiss;
            pdf_printf(left_margin, page_height,'L', "%li %s", options->newpops[pop], data->popnames[pop]);
            pdf_printf_ralign(col2, page_height,"1");
            pdf_printf_ralign(col3, page_height,"%li",numind);
            pdf_printf_right_next(left_margin+10, &page_height,"(%li)", nummiss);
        }
        else
        {
            nummiss = 0;
            numind =  data->numind[pop][0];
            //(options->randomsubset > 0 && options->randomsubset < data->numind[pop][0]) ? options->randomsubset : data->numind[pop][0];
            pdf_printf(left_margin, page_height, 'L', "%li %s",options->newpops[pop], data->popnames[pop]);
            pdf_printf_ralign(col2, page_height,"1");
            pdf_printf_ralign(col3, page_height,"%li", numind);
            pdf_advance(&page_height);
        }
        total[0] += numind;
        totalmiss[0] += nummiss;
        
        for(locus=1; locus< data->loci; locus++)
        {
            if (!strchr (SEQUENCETYPES, options->datatype) && options->datatype !='@')
            {
                nummiss = find_missing(data,pop,locus);
                numind = data->numalleles[pop][locus] - nummiss;
                pdf_printf_ralign(col2, page_height,"%li",locus+1);
                pdf_printf_ralign(col3, page_height,"%li",numind);
                pdf_printf_right_next(left_margin+10, &page_height,"(%li)", nummiss);
            }
            else
            {
                nummiss=0;
                numind = data->numind[pop][locus];
                pdf_printf_ralign(col2, page_height,"%li",locus+1);
                pdf_printf_ralign(col3, page_height,"%li",numind);
                pdf_advance(&page_height);
            }
            total[locus] += numind;
            totalmiss[locus] += nummiss;
        }
    }
    pdf_printf(left_margin, page_height, 'L',
               "Total of all populations");
    pdf_printf_ralign(col2, page_height,"1");
    pdf_printf_ralign(col3, page_height,"%li",total[0]);
    if (!strchr (SEQUENCETYPES, options->datatype))
    {
        pdf_printf_right_next(left_margin+10, &page_height,"(%li)", totalmiss[0]);
        for(locus=1; locus< data->loci; locus++)
        {
            pdf_printf_ralign(col2, page_height,"%li",locus+1);
            pdf_printf_ralign(col3, page_height,"%li",total[locus]);
            pdf_printf_right_next(left_margin+10, &page_height,"(%li)", totalmiss[locus]);
        }
    }
    else
    {
        pdf_advance(&page_height);
        for(locus=1; locus< data->loci; locus++)
        {
            pdf_printf_ralign(col2, page_height,"%li",locus+1);
            pdf_printf_ralign(col3, page_height,"%li",total[locus]);
            pdf_advance(&page_height);
        }
    }
    myfree(total);
    myfree(totalmiss);
    *orig_page_height = page_height;
    *orig_left_margin = left_margin;
}

///
// \param[in] *fmt format string identical to a fomrat string in printf()
/// \returns int returns the number of elements by counting the %
int count_elements(char *fmt)
{
    int element;
    int count = 0;
    for(element=0; element < (int) strlen(fmt); element++)
    {
        if(fmt[element]=='%')
            count++;
    }
    return count;
}


//
// Print  a line in a table in the PDF file at a specific height on the page and with a given width
// \param *page_height gets modified by this function, page_height specifies the Y coordinate
// \param *width contains a list of doubleing points numbers that specify the columns, negative numbers
// mean left-adjusted and positive numbers are right-adjusted
// \param *fmt format string identical to printf (number % needs to match the element in width)
// \param ...  parameters to print
void pdf_print_tableline(double *width, char *fmt, ...)
{
    boolean start=FALSE;
    boolean stop=FALSE;
    char fmtval;
    char this_fmt[LINESIZE];
    long ival;
    va_list ap;
    double dval;
    long fmti=0;
    char cval, *sval;
    double y = page_height;
    long ii=0;
    boolean left_aligned=FALSE;
    double offset = 55;
    va_start(ap, fmt);
    while (*fmt)
    {
        if(*fmt == '%')
        {
            if(width[ii] <= EPSILON)
            {
                left_aligned = TRUE;
                offset = 55 - width[ii++];
            }
            else
            {
                left_aligned = FALSE;
                offset = 55 + width[ii++];
            }
            
            start = TRUE;
            fmt++;
            this_fmt[0] = '%';
            fmti = 1;
            continue;
        }
        else
        {
            switch(*fmt)
            {
                case 'c':                       /* string */
                    stop = TRUE;
                    fmtval = 'c';
                    cval = (char) va_arg(ap, int);
                    this_fmt[fmti++] = fmtval;
                    this_fmt[fmti] = '\0';
                    if(left_aligned)
                        pdf_printf(offset, y, 'L', this_fmt, cval);
                    else
                        pdf_printf_ralign(offset, y,this_fmt, cval);
                    break;
                case 's':                       /* string */
                    stop = TRUE;
                    fmtval = 's';
                    sval = va_arg(ap, char *);
                    this_fmt[fmti++] = fmtval;
                    this_fmt[fmti] = '\0';
                    if(left_aligned)
                        pdf_printf(offset, y,'L', this_fmt, sval);
                    else
                        pdf_printf_ralign(offset, y,this_fmt, sval);
                    break;
                case 'i':                       /* int */
                    stop = TRUE;
                    fmtval = 'i';
                    ival = va_arg(ap, int);
                    this_fmt[fmti++] = fmtval;
                    this_fmt[fmti] = '\0';
                    if(left_aligned)
                        pdf_printf(offset, y,'L', this_fmt, ival);
                    else
                        pdf_printf_ralign(offset, y,this_fmt, ival);
                    break;
                case 'f':                       /* char */
                    stop = TRUE;
                    fmtval = 'f';
                    dval = va_arg(ap, double);
                    this_fmt[fmti++] = fmtval;
                    this_fmt[fmti] = '\0';
                    if(left_aligned)
                        pdf_printf(offset, y,'L', this_fmt, dval);
                    else
                        pdf_printf_ralign(offset, y,this_fmt, dval);
                    break;
                case '\n':
                    pdf_advance(&y);
                    break;
            }
            if(start)
            {
                if(!stop)
                {
                    if(strchr("1234567890.", *fmt)!=NULL)
                        this_fmt[fmti++] = *fmt;
                    if(*fmt == '-')
                        left_aligned = TRUE;
                }
                else
                {
                    start = FALSE;
                    stop = FALSE;
                }
            }
            fmt++;
        }
    }
    va_end(ap);
}


void pdf_print_allelelegend(double *column_width, long loci)
{
    long locus;
    char stemp[LINESIZE];
    sprintf(stemp, "%-s", (loci == 1 ? "locus" : "loci "));
    
    pdf_printf(column_width[0], page_height, 'L', "Indiv.");
    for(locus=1; locus < loci+1; locus++)
    {
        pdf_printf(column_width[locus], page_height, 'L', "%li",locus);
        if(column_width[locus+1]<column_width[locus])
        {
            pdf_advance(&page_height);
        }
    }
}

#define NOTFIRSTDATAPAGE (boolean) 0
#define FIRSTDATAPAGE (boolean) 1
///
/// print data header for allelic or sequence data
void pdf_print_dataheader (boolean first, char *title, double *column_width, long pop, world_fmt * world,
                      option_fmt * options, data_fmt * data)
{
    
    double w;
    double page_width;
    //double left_margin = 55;
    
    // setup new page and title
    if(first)
    {
        pdf_new_page("");
        HPDF_Page_SetFontAndSize (page, helvob, 18.0);
        w = (double) HPDF_Page_TextWidth(page, title);
        page_height = (double) HPDF_Page_GetHeight(page) - 75.0;
        page_width = (double) HPDF_Page_GetWidth(page);
        page_height -= 20.;
        pdf_print_contents_at((page_width - w)/2, page_height, title);
        HPDF_Page_SetRGBStroke(page, 0, 0, 0);
        page_height -= 20.;
        pdf_draw_line(50, page_height, page_width-50, page_height);
        page_height -= 20.;
        HPDF_Page_SetFontAndSize (page, helv, 10.0);
    }
    else
    {
        pdf_advance(&page_height);
        page_width = (double) HPDF_Page_GetWidth(page);
        pdf_draw_line(50, page_height, page_width-50, page_height);
        pdf_advance(&page_height);
    }
    pdf_printf_next(left_margin, &page_height, "\n%-s", data->popnames[pop]);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_advance(&page_height);
    if (!strchr (SEQUENCETYPES, options->datatype))
        pdf_print_allelelegend(column_width, world->loci);
    //    else
    //    pdf_print_seqlegend(left_margin, page_height, world->loci);
    
    pdf_advance(&page_height);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_advance(&page_height);
}



///
/// print the allele typ data in diploid format
// \param page_width page width
// \param page_height is manipulated and is always at the actual y coordinate on the page
// \param[in] margin is the margin that is used, either left or right margin
// \param[in] world contains all main data structures
// \param[in] data  contains all the biological data structures
// \param[in] options cintains all option structures
void pdf_print_alleledata (double margin, world_fmt * world, data_fmt * data, option_fmt * options)
{
    long top, ii;
    long pop, ind, locus;
    char stemp[LINESIZE];
    double w;
    double * column_width = (double *) mycalloc(world->loci+2, sizeof(double));
    double pixel;
    double total;
    double oldpage_height;
    double page_width = (double) HPDF_Page_GetWidth(page);
    
    // calculate column width for indvidual name and each locus
    // we use a cummulative sum to set the x-coord just right without further
    // calculations
    for (locus = 0; locus < data->loci; locus++)
    {
        for (pop = 0; pop < data->numpop; pop++)
        {
            if(options->randomsubset > 0 && options->randomsubset < data->numind[pop][locus])
            {
                top = options->randomsubset;
            }
            else
            {
                top = data->numind[pop][locus];
            }
            
            for (ii = 0; ii < top; ii++)
            {
                ind = data->shuffled[pop][locus][ii];
                
                w = (double)( 5. + (double) (strlen (data->yy[pop][ind][locus][0][0]) +
                                           strlen (data->yy[pop][ind][locus][1][0])));
                if(column_width[locus+2] < w)
                    column_width[locus+2] = w;
            }
        }
    }
    // calculate columnwidth for in page units for individual name and locus column
    column_width[1] = 5.0 + (double) (sprintf(stemp, "%-*.*s", (int) options->nmlength,
					    (int) options->nmlength, data->indnames[0][0][0]));
    w = (double) HPDF_Page_TextWidth(page, stemp);
    pixel = w / options->nmlength;
    column_width[0] = margin;
    column_width[1] = margin + w;
    total = margin + w;
    for (locus = 0; locus < data->loci; locus++)
    {
        column_width[locus+2] *= pixel;
        total += column_width[locus+2];
        if(total > (page_width - 55))
        {
            total = margin;
            column_width[locus+2] = total;
        }
        else
        {
            column_width[locus+2] = total;
        }
    }
    
    for (pop = 0; pop < data->numpop; pop++)
    {
        if(pop==0)
            pdf_print_dataheader (FIRSTDATAPAGE, "Allelic data", column_width, pop, world, options, data);
        else
            pdf_print_dataheader (NOTFIRSTDATAPAGE, "Allelic data", column_width, pop, world, options, data);
        
        for (ind = 0; ind < data->numind[pop][0]; ind++)
        {
            pdf_printf(column_width[0], page_height, 'L', "%-*.*s", (int) options->nmlength,
                       (int) options->nmlength, data->indnames[pop][ind][0]);
            for (locus = 0; locus < data->loci; locus++)
            {
                pdf_printf(column_width[locus+1], page_height, 'L', " %s.%-s",
                           data->yy[pop][ind][locus][0][0],
                           data->yy[pop][ind][locus][1][0]);
                if(column_width[locus+2]<column_width[locus+1])
                {
                    pdf_advance(&page_height);
                }
            }
            oldpage_height = page_height;
            pdf_advance(&page_height);
            if(oldpage_height < page_height)
                pdf_print_dataheader (NOTFIRSTDATAPAGE, "Allelic data", column_width, pop, world, options, data);
            
        }
    }
    myfree(column_width);
}




void
pdf_print_data (world_fmt * world, option_fmt * options, data_fmt * data)
{
    if (options->printdata)
    {
        switch (options->datatype)
        {
            case 'a':
            case 'b':
            case 'm':
                pdf_print_alleledata (left_margin, world, data, options);
                break;
            case 's':
            case 'n':
            case 'h':
            case 'u':
            case 'f':
                pdf_print_seqdata (left_margin, world, data, options);
                break;
            case '@':
                warning("PDF-printing of mixed data is not implemented ,yet\n");
                break;
        }
    }
}



void pdf_print_sequence(double right_margin, data_fmt *data, long locus, long pop, long ind)
{
    long site;
    double w = 0.;
    double wtot = 0.;
    char stemp[LINESIZE];
    //    w = 60.;
    w = 6.;
    for(site=0; site < data->seq[0]->sites[locus]; site+=1)
    {
        //        sprintf(stemp,"%-10.10s", data->yy[pop][ind][locus][0][site]);
        sprintf(stemp,"%1s", data->yy[pop][ind][locus][0][site]);
        HPDF_Page_SetFontAndSize (page, cour, 9.0);
        pdf_print_contents_at(left_margin + wtot, page_height, stemp);
        wtot += w;
        if((left_margin + wtot + w) > right_margin)
        {
            wtot = 0;
            pdf_advance(&page_height);
        }
    }
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
}

void pdf_print_seqdata (double margin, world_fmt * world, data_fmt * data, option_fmt * options)
{
    long *top;
    long ii;
    long pop, ind, locus=0;
    char stemp[LINESIZE];
    double w;
    double * column_width = (double *) mycalloc(world->loci+2, sizeof(double));
    //xcode   double pixel;
    double oldpage_height;
    double page_width = (double) HPDF_Page_GetWidth(page);
    top = (long *) mycalloc((data->loci * data->numpop),sizeof(long));
    // calculate columnwidth for in page units for individual name and locus column
    column_width[1] = 1.0 + (double) (sprintf(stemp, "%-*.*s", (int) options->nmlength,
					      (int) options->nmlength, data->indnames[0][0][0]));
    w = (double) HPDF_Page_TextWidth(page, stemp);
    //xcode pixel = w / options->nmlength;
    column_width[0] = margin;
    column_width[1] = margin + w;
    column_width[2] = margin + w + (double) 10;
    long maxtop = 0;
    for (pop = 0; pop < data->numpop; pop++)
    {
        for(locus=0;locus<data->loci;locus++)
        {
            top[pop*data->loci + locus] = max_shuffled_individuals(options, data, pop, locus);
            if(maxtop < top[locus])
                maxtop = top[locus];
        }
    }
    for (pop = 0; pop < data->numpop; pop++)
    {
        if(pop==0)
            pdf_print_dataheader (FIRSTDATAPAGE, "Sequence data", column_width, pop, world, options, data);
        else
            pdf_print_dataheader (NOTFIRSTDATAPAGE,  "Sequence data", column_width, pop, world, options, data);
        
        for (ii = 0; ii < maxtop; ii++)
        {
	  if(ii>=top[pop * data->loci])
	    {
	      continue;
	    }
	  else
	    {
	      ind = data->shuffled[pop][0][ii];
	      pdf_printf(column_width[0], page_height, 'L', "%-*.*s", (int) options->nmlength,
			 (int) options->nmlength, data->indnames[pop][ind][0]);
	    }
	  for (locus = 0; locus < data->loci; locus++)
            {
                if(ii>=top[pop * data->loci + locus])
                {
                    continue;
                }
                else
                {
                    ind = data->shuffled[pop][locus][ii];
                    pdf_printf(column_width[1], page_height, 'L', "%li",locus+1);
                    pdf_print_sequence( page_width - 55, data, locus, pop, ind);
                    pdf_advance(&page_height);
                }
            }
            oldpage_height = page_height;
            pdf_advance(&page_height);
            pdf_advance(&page_height);
            if(oldpage_height < page_height)
                pdf_print_dataheader (NOTFIRSTDATAPAGE,  "Sequence data", column_width, pop, world, options, data);
        }
    }
    myfree(column_width);
    myfree(top);
}

//void format_helper_vec(atl_fmt *atl, long pop, long locus, long maxrep, long numpop, long *fmt_vector)
//{
//  for (i=0; i < atl.params; i++)
//}

void format_helper(MYREAL * matrix, long len, int *fmt1, int *fmt2)
{
    long i;
    long count = 0;
    long maxcount = 0;
    MYREAL val;
    
    for(i=0;i<len; i++)
    {
        val = matrix[i];
        if(val>0.)
            count = (long)(log10(val)) + 1;
        if(count>maxcount)
            maxcount = count;
    }
    *fmt1 = 7;
    if(maxcount<2)
    {
        *fmt2 = 5;
    }
    else
    {
        *fmt2 = 2;
    }
}


///
/// print matrix line for various uses
/// - print linear shortened migration matrix in n x n matrix with - as diagonals
/// alowing for some pretty print using a format_helper() the columns are right-aligned
/// based on the width.
void pdf_print_matrix_line(long whichline, double width, long cols,long rows,MYREAL * matrix,boolean ismig)
{
    long i=0;
    long col, row = whichline;
    double offset;
    double page_width = (double) HPDF_Page_GetWidth(page);
    int fmt1=5, fmt2=3;
    long pad = ((ismig)?cols:0);
    format_helper(matrix+pad, cols * rows - pad ,&fmt1, &fmt2);
    offset = left_margin;
    if(ismig)
        i = cols;
    for(col = 0; col < cols; col++)
    {
        if(offset > page_width - 55)
        {
            offset = left_margin;
            pdf_advance(&page_height);
        }
        if(ismig && row==col)
        {
            pdf_printf_ralign(offset+width - 15, page_height, "-");
        }
        else
        {
            pdf_printf_ralign(offset + width - 5, page_height, "%*.*f", fmt1, fmt2, matrix[i]);
            i++;
        }
        offset += width;
    }
    pdf_advance(&page_height);
}

///
/// print text using a simple formatter that works on a single paramgraph
/// breaking it approx at the end of the line.
void pdf_print_comment(double lx, double *ly, char *this_text)
{
    pdf_printf(lx, *ly,'L', this_text);
    pdf_advance(ly);
}



/*
 * print header
 * ===========================================================================
 * === Titletext
 * ===========================================================================
 * === Population     Loci  Ln(L)   Theta    xNm [xNe mu] xx     xx     xx
 * xx     xx     xx xx.... -------------- ---- -------- --------
 * ----------------------------------------
 */
/// print the MLE table header
void pdf_print_result_header (double *lx, char *titletext, world_fmt * world)
{
    long p1, zz;
    double *ly = &page_height;
    double page_width;
    
    pdf_print_section_title(&page_width, ly, titletext);
    pdf_advance(&page_height);
    pdf_printf(lx[0], *ly,'L',"Population [x]");
    pdf_printf(lx[1], *ly,'L',"Loc.");
    pdf_printf(lx[2], *ly,'L',"Ln(L/L0)");
    symbol_Theta(lx[3], *ly,12,-1);
    //    pdf_printf(lx[3], *ly,'L',"Theta");
    if(world->numpop>1)
    {
        if(world->options->usem)
            pdf_printf(lx[4], *ly,'L',"M (m/mu) [+=receiving population");
        else
            pdf_printf(lx[4], *ly,'L',"xNm [+=receiving population");
    }
    pdf_advance(ly);
    pdf_printf(lx[3], *ly,'L',"[x Ne mu]");
    
    zz = 4;
    for (p1 = 0; p1 < world->numpop; p1++)
    {
        if (zz > 8)
        {
            zz = 4;
            pdf_advance(ly);
        }
        pdf_printf(lx[zz], *ly,'L',"%2li,+", p1 + 1);
        zz++;
    }
    pdf_advance(ly);
    pdf_draw_line(lx[0],*ly, page_width - 55, *ly);
    pdf_advance(ly);
}

/// print the string for a population
void pdf_print_popstring(double *lx, long pop, world_fmt *world, option_fmt *options, data_fmt *data)
{
  (void) world;
    double *ly = &page_height;
    //    char popstring[LINESIZE];
    if (options->readsum)
    {
        pdf_printf(lx[0], *ly, 'L', "%2li:",pop+1);
    }
    else
    {
        pdf_printf(lx[0],*ly, 'L', "%2li:%10.10s",pop+1, data->popnames[options->newpops[pop]-1]);
    }
}


/// print the replicate number
void pdf_print_replicate(double lx, world_fmt *world, long maxrep, long rep, long locus)
{
  (void) world;
    double *ly = &page_height;
    char repstring[LINESIZE];
    sprintf (repstring, "%li", rep + 1);
    pdf_printf(lx, *ly, 'L', "%li %s", locus + 1, maxrep > 1 ? (rep == maxrep - 1 ? " A" : repstring) : "  ");
}

/* ML material needs to be excised
/// print the MLE table content for each population
void pdf_print_result_population (double *lx, long pop, world_fmt * world,
                             option_fmt * options, data_fmt * data)
{
    double *ly = &page_height;
    long skipped = 0, locus;
    long maxrep = world->options->replicate ?
    (world->options->replicatenum > 0 ?
     world->options->replicatenum + 1 : world->options->lchains + 1) : 1;
    long rep;
    pdf_print_popstring(lx, pop, world, options, data);
    for (locus = 0; locus < world->loci; locus++)
    {
        if (world->data->skiploci[locus])
        {
            skipped++;
            continue;
        }
        for (rep = 0; rep < maxrep; rep++)
        {
            pdf_print_replicate(lx[1], world, maxrep, rep, locus);
            pdf_printf(lx[2], *ly, 'L', "% 8.3f ",
                       world->atl[rep][locus].param_like);
            pdf_print_result_param (lx,  world->atl[rep][locus].param,
                                    world->numpop, pop, world->options->usem);
        }
    }
    if (world->loci - skipped > 1)
    {
        pdf_printf(lx[1], page_height, 'L', "All ");
        //locus is exactly world->loci
        // re is always one because we have only replication of single locus chains
        pdf_printf(lx[2],page_height, 'L', "% 8.3f ", world->atl[0][locus].param_like);
        pdf_print_result_param (lx,  world->atl[0][locus].param,
                                world->numpop, pop, world->options->usem);
    }
    // FPRINTF(world->outfile,"%s\n",sline);     
}


/// print the parameter for the MLE parameter printout
void
pdf_print_result_param (double *lx,  MYREAL *param, long numpop, long pop,
                        boolean usem)
{
    char    temp[LINESIZE];
    long    i;
    long    zz;
    long    msta = mstart (pop, numpop);
    long    msto = mend (pop, numpop);
    double  *ly   = &page_height;
    //int   fmtint = 8;
    //int fmtdouble = 5;
    //int digits;
    MYREAL  tmp  = 0;
    
    // population size
    if (param[pop] <= SICK_VALUE)
        pdf_printf (lx[3], *ly, 'L', "-");
    else
    {
        nice_element(param[pop], temp, 0.0001, 10, 100, 4, 2, '\0');
        pdf_printf (lx[3], *ly, 'L', "%s", temp);
        //        if (param[pop] < 0.0001 && param[pop] > 0)
        // pdf_printf (lx[3], *ly, 'L', "%3.2e ", param[pop]);
        // else
        // {
        // digits = (int) log10(param[pop]);
        // if(digits>3)
        // {
        // fmtdouble=0;
        // fmtint=digits;
        // }
        // pdf_printf (lx[3], *ly, 'L', "%*.*f ", fmtint, fmtdouble, param[pop]);
        // }
         ///
    }
    
    // migration rate
    zz=4;
    for (i = msta; i < msto; i++)
    {
        if (zz > 8)
        {
            zz = 4;
            pdf_advance(ly);
        }
        if (pop == i - msta)
        {
            pdf_printf (lx[zz], *ly, 'L', "-");
            zz++;;
        }
        if ((param[i] <= SICK_VALUE) || (param[pop] <= SICK_VALUE))
            pdf_printf (lx[zz], *ly, 'L', "-");
        else
        {
            if (usem)
            {
                //tmp = param[i];
                nice_element(param[i], temp, 0.001, 100, 1000, 3, 2, '\0');
                //if (tmp < 0.0001)
                pdf_printf (lx[zz], *ly, 'L', "%s",temp);
                //else
                //pdf_printf (lx[zz], *ly, 'L', "%7.4f ", tmp);
            }
            else
            {
                tmp = param[pop] * param[i];
                nice_element(tmp, temp, 0.0001, 10, 100, 4, 2, '\0');
                pdf_printf (lx[zz], *ly, 'L', "%s",temp);
                //		if (tmp < 0.00001)
                // pdf_printf (lx[zz], *ly, 'L', " 0.0000 ");
                // else
                // pdf_printf (lx[zz], *ly, 'L', "%7.4f ", tmp);
                 //
            }
        }
        zz++;
    }
    if (pop == numpop - 1)
        pdf_printf (lx[i-msta+4], *ly, 'L',"-");
    pdf_advance(ly);
}


void
pdf_print_results (world_fmt ** universe, option_fmt * options, data_fmt * data)
{
    double lx[]={55,140,170,220,270,320,370,420,470,520};
    
    long pop;
    double *ly = &page_height;
    //xcode    FILE *outfile;
    world_fmt *world = universe[0];
    worldoption_fmt *wopt = world->options;
    char sch[10], lch[10], cva[50];
    long rep = world->loci > 1 ? 0 : (wopt->replicate ? world->repstop : 0);
    //xcode outfile = world->outfile;
    if (options->schains == 1)
        strcpy (sch, "chain");
    else
        strcpy (sch, "chains");
    if (options->lchains == 1)
        strcpy (lch, "chain");
    else
        strcpy (lch, "chains");
    pdf_print_result_header (lx, "Maximum Likelihood estimates", world);
    for (pop = 0; pop < world->numpop; pop++)
    {
        pdf_print_result_population (lx, pop, world, options, data);
    }
    pdf_advance(ly);
    pdf_printf(lx[0],*ly,'L', "Comments:");
    pdf_advance(ly);
    pdf_printf(lx[0],*ly,'L', "The x is 1, 2, or 4 for mtDNA, haploid, or diploid data, respectively");
    pdf_advance(ly);
    pdf_printf(lx[0],*ly,'L', "There were %li short %s (%li used trees out of sampled %li)",
               options->schains, sch, options->ssteps, options->sincrement * options->ssteps);
    pdf_advance(ly);
    pdf_printf(lx[0],*ly, 'L', "and %li long %s (%li used trees out of sampled %li)\n",
               options->lchains, lch, options->lsteps,
               options->lincrement * options->lsteps);
    pdf_advance(ly);
    if (wopt->heating)
    {
        if(options->adaptiveheat!=NOTADAPTIVE)
        {
            if(options->adaptiveheat==STANDARD)
            {
                pdf_printf(lx[0],*ly,'L',
                           "Adaptive heating with %li chains was active",options->heated_chains);
            }
            else
            {
                pdf_printf(lx[0],*ly,'L',
                           "Bounded adaptive heating with %li chains was active",options->heated_chains);
            }
            pdf_advance(ly);
            pdf_printf(lx[0],*ly,'L',"check Log file (if present) for average temperatures");
            
            pdf_printf(lx[0],*ly,'L',"Average last chain temp: 1.0",
                       options->heated_chains);
            for(pop=1;pop<options->heated_chains;pop++)
                pdf_printf(lx[pop+2],*ly,'L', ", %f", universe[pop]->averageheat);
            pdf_advance(ly);
        }
        else
        {
            pdf_printf(lx[0],*ly,'L', "Static heating with %li chains was active\n",options->heated_chains);
            pdf_advance(ly);
        }
    }
    if (options->gamma)
    {
        if (world->atl[rep][world->loci].param[world->numpop2] < 10e-9)
            strcpy (cva, "0");
        else
            sprintf (cva, "%f",
                     sqrt (1. / world->atl[rep][world->loci].param[world->numpop2]));
        pdf_printf(lx[0],*ly,'L',"With shape parameter Alpha=%g ([1/CV(mu)]^2; CV(mu)=%s)",
                   world->atl[rep][world->loci].param[world->numpop2],
                   cva);
        pdf_advance(ly);
    }
    if (world->options->replicate)
    {
        if (world->repkind == MULTIPLECHAIN)
            pdf_printf(lx[0],*ly, 'L', "COMBINATION OF ALL LONG CHAINS");
        else
            pdf_printf(lx[0],*ly, 'L', "COMBINATION OF %li MULTIPLE RUNS",
                       world->options->replicatenum);
        pdf_advance(ly);
    }
    if (world->atl[rep][world->loci].normd > LOCI_NORM)
    {
        pdf_printf(lx[0],*ly,'L',"[Last maximization needed %li cycles of maximal %i,",
                   world->atl[rep][world->loci].trials,
                   NTRIALS);
        pdf_advance(ly);
        pdf_printf(lx[0],*ly,'L',"Norm(first derivatives)=%f (Normal stopping criteria is < %f)]",
                   world->atl[rep][world->loci].normd, LOCI_NORM);
    }
    pdf_advance(ly);
    pdf_advance(ly);
    pdf_print_citation("Likelihood inference", world);
}
*/

void    pdf_print_correlation_table (world_fmt * world, option_fmt * options, data_fmt * data)
{
    const long nn = world->numpop2 + world->bayes->mu;
    const double offset = 60;
    long i;
    long j;
    double right_margin = (double) HPDF_Page_GetWidth(page) - 2.0*55.0 - offset;
    //double left_margin = 55;
    double lx=left_margin + 80;
    long frompop, topop;
    long locus;
    long loci = world->loci > 1 ? world->loci+1 : world->loci;
    pdf_advance(&page_height);
    pdf_printf_next(left_margin, &page_height, "Correlation tables:");
    pdf_advance(&page_height);
    for (locus=0; locus < loci; locus++)
    {
        pdf_printf (left_margin, page_height, 'L', "Locus %li",locus+1);
        for (i = 0; i < nn; i++)
        {
            if(i<world->numpop)
                symbol_Theta(left_margin, page_height, 12, i+1);
            else
                if (i==world->numpop2)
                    symbol_R(left_margin, page_height, 12, -1);
                else
                {
                    m2mm(i,world->numpop,&frompop,&topop);
                    symbol_M(left_margin, page_height, 12, frompop+1, topop+1, world->options->usem);
                }
            
            lx += offset;
            if(lx > right_margin)
            {
                lx = left_margin + 80;
                pdf_advance(&page_height);
            }
        }
    }
    for (i = 0; i < data->numpop; i++)
    {
        lx = left_margin + 80;
        pdf_advance(&page_height);
        pdf_printf (left_margin, page_height, 'L', "%3li %-15.15s",options->newpops[i],  data->popnames[i]);
        for (j = 0; j < data->numpop; j++)
        {
            pdf_printf(lx, page_height, 'L', " %10.4f ", data->ogeo[j][i]);
            lx += offset;
            if(lx > right_margin)
            {
                lx = left_margin + 80;
                pdf_advance(&page_height);
            }
        }
    }
    pdf_advance(&page_height);
    pdf_advance(&page_height);
}


///
/// pretty print for tables
/// takes a lower and upper bound below and above it uses scientific notation, the range in
/// between prints fixed point notation below mid with low_mid_digits and above mid with mid_upper_digits.
long nice_element(MYREAL param, char *element, MYREAL lower, MYREAL mid, MYREAL upper,
                  int low_mid_digits, int mid_upper_digits, char delimiter)
{
    long position = 0;
    if (param < lower)
        position = sprintf(element,"%.2e", param);
    else
    {
        if (param > upper)
            position = sprintf(element,"%.2e", param);
        else
        {
            if(param > mid)
                position = sprintf(element,"%.*f", mid_upper_digits,param);
            else
                position = sprintf(element,"%.*f", low_mid_digits,param);
        }
    }
    if(delimiter != '\0')
        position += sprintf(element + position,"%c",delimiter);
    return position;
}



void method_set(double lx, char method)
{
    switch (method)
    {
        case 'p':
            pdf_printf(lx, page_height, 'L', "Parameters are evaluated at percentiles using bisection method (slow, but exact).");
            break;
        case 's':
            pdf_printf(lx, page_height, 'L',"Parameters are evaluated at percentiles");
            pdf_advance(&page_height);
            pdf_printf(lx, page_height, 'L',"using cubic splines of profiled parameter (faster, but not so exact).");
            break;
        case 'd':
            pdf_printf(lx, page_height, 'L',"Parameters are evaluated at pre-defined values\n");
            break;
        case 'q':
            pdf_printf(lx, page_height, 'L', "Parameters are evaluated assuming complete independence\n");
            break;
        case 'f':
            pdf_printf(lx, page_height, 'L', "Parameters are evaluated assuming complete independence");
            pdf_advance(&page_height);
            pdf_printf(lx, page_height, 'L',"and then once maximized at the found profiled parameter value");
            break;
    }
    pdf_advance(&page_height);
}

void pdf_table_footnote(double lx, long failed, boolean percentiles)
{
    char * foot1 = "If the percentile column is marked with **** then the convergence to the percentile value failed";
    char * foot2 = "The values are still CORRECT but not at the percentile of the profile parameter";
    char * foot3 = "Values with a * are NOT at the percentiles!";
    char * foot4 = "The convergence to the percentile value FAILED.";
    if(failed>0)
    {
        pdf_printf(lx, page_height, 'L', percentiles ? foot3 : foot1);
        pdf_advance(&page_height);
        pdf_printf(lx, page_height, 'L', percentiles ? foot4 : foot2);
        pdf_advance(&page_height);
    }
}

///
/// translate buffer table (syntax similar to LaTeX) into table header and table elements
void  translate_buffer_table(long cols, long rows, char **thebuffer, char **header, char ***elements)
{
    char *temp;
    char *buffer;
    char *savebuffer;
    long z = 0;
    long r = 0;
    //  long position=0;
    
    buffer = (char *) mycalloc(strlen(*thebuffer)+1,sizeof(char));
    savebuffer= buffer;
    strcpy(buffer,*thebuffer);
    while((temp=strsep(&buffer,"&"))!=NULL && buffer!=NULL )
    {
        if(temp[0]=='%')
            break;
        else
        {
            sprintf(header[z++],"%s",temp);
        }
    }
    z=0;
    r=0;
    while((temp=strsep(&buffer,"&"))!=NULL && buffer!=NULL)
    {
        if(r==rows-1 && z>=cols)
        {
            //	  position = (long) strlen(strstr(thebuffer,buffer))+1;
            // memcpy(*thebuffer,buffer,(strlen(buffer)+1)*sizeof(char));
            break;
        }
        if(temp[0]=='%')
            break;
        if(temp[0]=='@')
        {
            z=0;
            r++;
        }
        else
        {
            if(z<cols)
            {
                sprintf(elements[r][z],"%s",temp);
                z++;
            }
            else
            {
                if(r<rows)
                    warning("column counter exceeded available columns: z=%li temp=%s\n",z, temp);
                else
                    break;
            }
        }
    }
    myfree(savebuffer);
    //  return position;
}
///
/// extract column out of  buffer table (syntax similar to LaTeX) into table header and table elements
void  extract_column_buffer_table(long col, long cols, long rows, char **thebuffer, double *x)
{
    char *temp;
    char *buffer;
    char *savebuffer;
    long z = 0;
    long r = 0;
    //  long position=0;
    
    buffer = (char *) mycalloc(strlen(*thebuffer)+1,sizeof(char));
    savebuffer= buffer;
    strcpy(buffer,*thebuffer);
    while((temp=strsep(&buffer,"&"))!=NULL && buffer!=NULL )
    {
        if(temp[0]=='%')
            break;
        else
        {
            z++;
            //            sprintf(header[z++],"%s",temp);
        }
    }
    z=0;
    r=0;
    while((temp=strsep(&buffer,"&"))!=NULL && buffer!=NULL)
    {
        if(r==rows-1 && z>=cols)
        {
            //	  position = (long) strlen(strstr(thebuffer,buffer))+1;
            // memcpy(*thebuffer,buffer,(strlen(buffer)+1)*sizeof(char));
            break;
        }
        if(temp[0]=='%')
            break;
        if(temp[0]=='@')
        {
            z=0;
            r++;
        }
        else
        {
            if(z==col)
	      x[r] = (double) atof(temp);
            
            if(z<cols)
            {
                //sprintf(elements[r][z],"%s",temp);
                z++;
            }
            else
            {
                if(r<rows)
                    warning("column counter exceeded available columns: z=%li temp=%s\n",z, temp);
                else
                    break;
            }
        }
    }
    myfree(savebuffer);
    //  return position;
}


// generic table generator package
// containing functions:
// pdf_table()
// find_col_width()
// define_col_start()
// pdf_print_table_header()

///
/// finds the width of all columns given the elements of the table
void  find_col_width(int cols, int rows, char ***elements, char **header, double *col_widths)
{
    double w;
    double keepw;
    int row;
    int col;
    
    for(col=0;col < cols; col++)
    {
        keepw = (double) HPDF_Page_TextWidth(page, header[col]);
        for(row=0;row < rows; row++)
        {
            w = (double) HPDF_Page_TextWidth(page, elements[row][col]);
            if(w > keepw)
                keepw = w;
        }
        col_widths[col] = keepw;
    }
}

///
/// align column
double align_column(char position, double cw, double col_leftmargin)
{
    double loc = 0.;
    switch(position)
    {
        case 'R':
            loc = col_leftmargin + cw;
            break;
        case 'C':
            loc = col_leftmargin + cw/2.0;
            break;
        case 'L':
        default:
            loc = col_leftmargin;
            break;
    }
    return loc;
}

///
/// sets X-coordinates to start the columns
void  define_col_start(double cols, double * col_widths, int col_overflow, char *position, double page_width, double separator, double *col_starts)
{
    int col;
    char pos;
    if(position==NULL)
        pos='L';
    else
        pos=position[0];
    col_starts[0] = align_column(pos, col_widths[0], left_margin);
    for(col=1; col < cols; col++)
    {
        if(position==NULL)
            pos='L';
        else
            pos=position[col];
        col_starts[col] =  align_column(pos, col_widths[col], separator + col_starts[col-1] + col_widths[col-1]);
        if((col_starts[col] + col_widths[col]) > page_width) //correct for 'L', but questionable for 'R'
        {
            if(col_overflow==0)
                col_starts[col] = col_starts[0];
            else
            {
                if(col_overflow <= col)
                    col_starts[col] = col_starts[col_overflow];
                else
                    col_starts[col] = col_starts[0];
            }
            //	  fprintf(stdout,">>> %f %f %f\n",col_starts[col],col_widths[col], separator);
        }
    }
}

///
/// prints the table header at column positions
void  pdf_print_table_header(char * position, int cols, double *col_starts,char **header)
{
    int col;
    char pos;
    long oldcolstarts = -1;
    for(col=0; col < cols; col++)
    {
        if(position==NULL)
            pos='L';
        else
            pos=position[col];
        if(col_starts[col] < oldcolstarts)
        {
            pdf_advance(&page_height);
        }
        pdf_print_symbol(col_starts[col], 10, pos, header[col]);
        oldcolstarts = (long) col_starts[col];
    }
}

void pdf_print_symbol_no(char value, char* symbolstring, double lx, char pos)
{
  (void) value;
    //if (value!='?')
    pdf_printf(lx, page_height, pos, "%s", symbolstring);
    //    else
    //pdf_printf(lx, *page_height, pos, "%s", "?");
}


void pdf_print_symbol(double lx, int fontsize, char pos, char *symbolstring)
{
    int topop;
    int frompop;
    char value = '?';
    if (symbolstring != NULL)
    {
        value = symbolstring[0];
    }
    switch(value)
    {
        case 'T':
            if(symbolstring[1]=='h')
            {
                sscanf(symbolstring,"Theta_%i",&topop);
                symbol_Theta(lx, page_height, fontsize, topop);
            }
            else
                pdf_print_symbol_no(value, symbolstring, lx, pos);
            break;
        case 'Q':
            if(symbolstring[1]=='_')
            {
                sscanf(symbolstring,"Q_%i",&topop);
                symbol_Theta(lx, page_height, fontsize, topop);
            }
            else
                pdf_print_symbol_no(value, symbolstring, lx,pos);
            break;
        case 'M':
            if(symbolstring[1]!='_')
                pdf_printf(lx, page_height, pos, "%s", symbolstring);
            else
            {
                sscanf(symbolstring,"M_%i->%i",&frompop,&topop);
                symbol_M(lx, page_height, fontsize, frompop, topop, TRUE);
            }
            break;
        case 'x':
            if(symbolstring[1]=='N')
            {
                sscanf(symbolstring,"xNm_%i->%i",&frompop,&topop);
                symbol_M(lx, page_height, fontsize, frompop, topop, TRUE);
            }
            else
                pdf_print_symbol_no(value, symbolstring, lx, pos);
            break;
        case '_':
            switch (symbolstring[1])
        {
            case 'T':
                pdf_printf(lx, page_height, pos, "%s", symbolstring+1);//should print TOTAL
                break;
            case 'H':
                symbol_Hexp(lx, page_height, fontsize);
                break;
        }
            break;
        default:
            pdf_print_symbol_no(value, symbolstring, lx,pos);
    }
}

//
// generic table generator using header and elements arrays
// \param double left_margin left edge of table
// \param double * page_height page height
// \param int cols number of columns in table if there are too many columns then new line
// \param int rows number of rows in table, if a new page is needed then header is repeated
// \param char *** elements all table elements, currently no formatting of these
// \param char ** header   header row
// \param char ** header2  second header row NULL of not needed
// \param char *position position of all elements: L=left, R=right, C=center
// \param int col_overflow when there are too many columns this is the column to restart
long pdf_table2(int cols, int rows, char **header, char **header2, char ***elements, char *position, int col_overflow, double separator)
{
    boolean new_page=FALSE;
    
    int row;
    int col;
    long location=0;
    
    char pos;
    
    double *col_widths;
    double page_width = (double) HPDF_Page_GetWidth(page);
    double right_margin = page_width - 55;
    double *col_starts;
    double oldcol = DBL_MAX - 1.0;
    
    col_widths = (double *) mycalloc(cols,sizeof(double));
    col_starts = (double *) mycalloc(cols,sizeof(double));
    
    find_col_width(cols, rows, elements, header, col_widths);
    define_col_start(cols, col_widths, col_overflow, position, right_margin, separator, col_starts);
    
    pdf_print_table_header(position, cols, col_starts,header);
    if (header2!=NULL)
    {
        pdf_advance(&page_height);
        pdf_print_table_header(position, cols, col_starts,header2);
    }
    new_page = pdf_advance(&page_height);
    pdf_draw_line(55,page_height+5.0,right_margin, page_height+5.0);
    
    for(row=0; row< rows; row++)
    {
        if(new_page)
        {
	  pdf_print_table_header(position, cols,col_starts,header);
            if (header2!=NULL)
            {
                pdf_advance(&page_height);
                pdf_print_table_header(position, cols, col_starts,header2);
            }
            pdf_advance(&page_height);
            pdf_draw_line(55.,page_height,right_margin, page_height);
            new_page = pdf_advance(&page_height);
        }
        for(col=0;col < cols; col++)
        {
            if(col_starts[col] < oldcol)
            {
                new_page = pdf_advance(&page_height);
                if(new_page)
                {
                    pdf_print_table_header(position, cols,col_starts,header);
                    if (header2!=NULL)
                    {
                        pdf_advance(&page_height);
                        pdf_print_table_header(position, cols, col_starts,header2);
                    }
                    pdf_advance(&page_height);
                    pdf_draw_line(55.,page_height,right_margin, page_height);
                    new_page = pdf_advance(&page_height);
                }
            }
            if(position==NULL)
                pos='L';
            else
                pos=position[col];
            pdf_print_symbol(col_starts[col],10,pos,elements[row][col]);
            //            pdf_printf(col_starts[col],*page_height, pos, "%s",elements[row][col]);
            oldcol = col_starts[col];
        }
    }
    pdf_advance(&page_height);
    myfree(col_widths);
    myfree(col_starts);
    return location;
}

//
// generic table generator using a buffer that contains header and elements, these
// are then split into header and elements, this approach allows to use the existing
// machinery for MPI.
// \param double left_margin left edge of table
// \param double * page_height page height
// \param int cols number of columns in table if there are too many columns then new line
// \param int rows number of rows in table, if a new page is needed then header is repeated
// \param char *** elements all table elements, currently no formatting of these
// \param char ** header   header row
// \param char *position position of all elements: L=left, R=right, C=center
// \param int col_overflow when there are too many columns this is the column to restart
void pdf_table(int cols, int rows, char **buffer, char *position, int col_overflow, double separator)
{
    int row;
    //  long location=0;
    char **header;
    char ***elements;
    
    charvec2d(&header, cols, LINESIZE);
    elements = (char ***) mycalloc(rows, sizeof(char **));
    for(row=0; row < rows; row++)
    {
        charvec2d(&(elements[row]),cols, LINESIZE);
    }
    
    translate_buffer_table(cols, rows, buffer, header, elements);
    
    pdf_table2(cols, rows, header, NULL, elements, position, col_overflow, separator);
    
    free_charvec2d(header);
    for(row=0; row < rows; row++)
        free_charvec2d(elements[row]);
    myfree(elements);
    // return location;
}




///
/// plot event histograms
void
pdf_event_histogram(long loci, long numparams,  world_fmt *world)
{
    double   page_width;
    double   lx;
    double   ly;
    double   ph;
    double   delta;
    double   w;
    double * binning;
    
    long    eventbinnum_allmax = 0L;
    long    loc;
    long    i;
    long    i0;
    long    j;
    long    bins;
    long    frompop;
    long    topop;
    long  * eventbinnum;
    
    char  * set50;
    char  * set95;
    char    title[LINESIZE];
    //double total;
    
    double ** eventbins_all;
    
    duo ** eventbins;
    
    
    if (loci > 1)
        sprintf(title,"Summary of events through time over all loci");
    else
        sprintf(title,"Events through time");
    
    // add a new page so that we can print at least four histograms
    pdf_new_page("");
    HPDF_Page_SetFontAndSize (page, helvob, 16.0);
    w = (double) HPDF_Page_TextWidth(page, title);
    
    /* Start to print text. */
    page_height = (double) HPDF_Page_GetHeight(page);
    page_width = (double) HPDF_Page_GetWidth(page);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    HPDF_Page_SetRGBStroke(page, 0, 0, 0);
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
    // first histogram position
    page_height -= 150;
    lx = 100;
    ly = page_height;
    
    //    for (loc = 0; loc < loci; loc++)
    //  {
    loc = loci > 1 ? loci : 0;
    eventbinnum = world->mighistloci[loc].migeventbinnum;
    for(i0=0; i0< numparams; i0++)
    {
      if(shortcut(i0,world,&i))
	{
	  continue;
	}
      
        if(eventbinnum[i] > eventbinnum_allmax)
            eventbinnum_allmax = eventbinnum[i];
    }

    doublevec2d(&eventbins_all,numparams,eventbinnum_allmax);
    //total = 0.0;
    for (loc = 0; loc < loci; loc++)
    {
        eventbinnum = world->mighistloci[loc].migeventbinnum;
        eventbins = world->mighistloci[loc].migeventbins;
        for(i0=0; i0< numparams; i0++)
        {
	  if(shortcut(i0,world,&i))
	    {
	      continue;
	    }
	  for(j=0; j < eventbinnum[i]; j++)
	    {
	      if(eventbins[i][j][0] > 0.0f)
		{
		  eventbins_all[i][j] += (double) eventbins[i][j][0];
		  //total += (MYREAL) eventbins[i][j][0];
		}
	    }
        }
    }
    bins = eventbinnum_allmax;
    
    set50 = (char *) mycalloc(bins+1, sizeof(char));
    set95 = (char *) mycalloc(bins+1, sizeof(char));
    memset (set50, '0', (size_t) bins * sizeof(char));
    memset (set95, '0', (size_t) bins * sizeof(char));
    
    binning = (double *) mycalloc (bins+1, sizeof (double));
    
    delta = (double) world->mighistloci[0].eventbinsize;
    binning[0] = 0.5 * delta;
    for (i = 1; i < bins; i++)
        binning[i] = delta + binning[i - 1];
    
    for (i0 = (world->options->mighist_all ? 0 : world->numpop); i0 < numparams; i0++)
    {
        if(shortcut(i0,world,&i))
	{
	  continue;
	}
	if(i<world->numpop2)
	  m2mm(i, world->numpop, &frompop, &topop);
	else
	  {
	    if (i0==world->numpop2 && world->bayes->mu)
	      {
		warning("histogram of expected rate changes not implemented");
		continue;
	      }
	    else
	      d2mm(i, world, &frompop, &topop);
	  }
        if(frompop==topop)
        {
            pdf_print_contents_at(lx-30,ly+125,"Freq. for ");
            symbol_Theta(lx+12, ly+125,12,frompop+1);
        }
        else
        {
            pdf_print_contents_at(lx-30,ly+125,"Freq. for ");
	    if (i0<world->numpop2)
	      symbol_M(lx+12, ly+125, 12, frompop+1, topop+1, world->options->usem);
	    else
	      {
		if (i0==world->numpop2 && world->bayes->mu)
		  symbol_R(lx+12, ly+125, 12,-1);
		else
		  symbol_D(lx+12, ly+125, 12, frompop+1, topop+1);
	      }
        }
        pdf_print_contents_at(lx+90,ly-25, "Time [scaled by mutation rate / site / generation]");
        HPDF_Page_SetFontAndSize (page, helv, 10.0);
        pdf_histogram(eventbins_all[i], set50, set95, bins, delta,
                      0., -999, lx, ly, page_width - 55 - lx, 116,TRUE,NULL);
        page_height -= 160;
        if(i < (numparams-1))
	  {
            pdf_page_advance_or_not(&page_height, 50);
            ph = (double) HPDF_Page_GetHeight(page) - 55.0 - LINESTRETCH;
            if(page_height >= ph)
	      page_height -= 150;
            ly = page_height;
	  }
    }
    myfree(binning);
    myfree(set50);
    myfree(set95);
    free_doublevec2d(eventbins_all);
}


///
/// plot skyline histograms
void
pdf_skyline_histogram(long loci, long numparams,  world_fmt *world, boolean enlarged)
{
  boolean   mu = world->bayes->mu;
    double   page_width;
    double   lx;
    double   ly;
    double   ph;
    double   delta;
    double   w;
    double * binning;
    double **confidence;
    double   c;
    double   lasttimebin;
    double   upperlimit;
    double maxtime;
    MYREAL *eventbin1max;
    //    double  confidencesum;
    long    clong;
    long    eventbinnum_allmax = 0L;
    long    loc;
    long    i;
    long    i0;
    long    j;
    long    bins;
    long    adjustbins;
    long    pop;
    long    frompop;
    long    topop;
    long  * eventbinnum;
    long    np2   = world->numpop2;
    long    nspec = world->species_model_size;
    long    npall = np2 + world->bayes->mu + nspec * 2;

    char  * set50;
    char  * set95;
    char    title[LINESIZE];
    
    
    double ** eventbins_all;
    MYREAL **std;
    
    tetra ** eventbins;
    species_fmt *s;
    
    if (loci > 1)
    {
        sprintf(title,"Summary of parameter values through %s over all loci",
                enlarged ? "RECENT time" : "time");
    }
    else
    {
        sprintf(title,"Parameter values through %s",enlarged ? "RECENT time" : "time");
    }
    // add a new page so that we can print at least four histograms
    pdf_new_page("");
    HPDF_Page_SetFontAndSize (page, helvob, 16.0);
    w = (double) HPDF_Page_TextWidth(page, title);
    
    /* Start to print text. */
    page_height = (double) HPDF_Page_GetHeight(page);
    page_width = (double) HPDF_Page_GetWidth(page);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    HPDF_Page_SetRGBStroke(page, 0, 0, 0);
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
    // first histogram position
    page_height -= 150;
    lx = 100;
    ly = page_height;
    
    //    for (loc = 0; loc < loci; loc++)
    // {
    loc = (loci > 1) ? loci : 0;
    eventbinnum = world->mighistloci[loc].eventbinnum;
    for(i0=0; i0< npall; i0++)
    {
        if(shortcut(i0,world,&i))
	{
	  continue;
	}
      
        if(eventbinnum[i] > eventbinnum_allmax)
            eventbinnum_allmax = eventbinnum[i];
    }
    // }
    doublevec2d(&eventbins_all,npall, eventbinnum_allmax+1);
    doublevec2d(&std,npall, eventbinnum_allmax+1);
    doublevec2d(&confidence,npall, eventbinnum_allmax+1);
    eventbin1max = (MYREAL *) mycalloc(npall,sizeof(MYREAL));
    loc = (loci > 1) ? loci : 0;
    eventbinnum = world->mighistloci[loc].eventbinnum;
    eventbins = world->mighistloci[loc].eventbins;
    for(i0=0; i0< numparams; i0++)
    {
        if(shortcut(i0,world,&i))
	{
	  continue;
	}    
        for(j=0; j < eventbinnum[i]; j++)
        {
            if(eventbins[i][j][1] > 0.0f)
            {
	      eventbins_all[i][j] = (double) eventbins[i][j][0];
	      std[i][j] = (MYREAL) eventbins[i][j][2] * 1.96; // standard error * 1.96
	      //old e...[2] changed to stderr: std[i][j] = (MYREAL) eventbins[i][j][2]; // standard deviation
	      if((double) eventbins[i][j][1] > eventbin1max[i])
		eventbin1max[i] = (double) eventbins[i][j][1];
            }
        }
    }
    
    for(i0=0; i0< numparams; i0++)
      {
	if(shortcut(i0,world,&i))
	  {
	    continue;
	  }    
        for(j=0; j < eventbinnum[i]; j++)
        {
            //c = MAX(0.0,MIN(1.0,(double) (eventbins[i][j][1])));
            clong = (long) eventbins[i][j][4];
            if(clong < 100 )
                c = 0.01;
            else if(clong < 400)
                c = 0.1;
            else if (clong < 1600)
                c = 0.5;
            else if (clong < 3200)
                c = 0.70;
            else if (clong < 6400)
                c = 0.8;
            else if (clong < 12800)
                c = 0.9;
            else
                c = 1.0;
            
            confidence[i][j] = 1.0 - c;
        }
    }
    
    bins = eventbinnum_allmax;
    
    set50 = (char *) mycalloc(bins+1, sizeof(char));
    set95 = (char *) mycalloc(bins+1, sizeof(char));
    memset (set50, '0', (size_t) bins * sizeof(char));
    memset (set95, '0', (size_t) bins * sizeof(char));
    
    binning = (double *) mycalloc (bins+1, sizeof (double));
    
    delta = (double) world->mighistloci[0].eventbinsize;
    binning[0] = 0.5 * delta;
    for (i = 1; i < bins; i++)
        binning[i] = delta + binning[i - 1];
    
    for (i0 = (world->options->mighist_all ? 0 : world->numpop); i0 < numparams; i0++)
    {
      frompop = -1;
      topop = -1;
      
      if(shortcut(i0,world,&i))
	{
	  continue;
	}

	if(i<world->numpop2)    
	  {
	    m2mm(i, world->numpop, &frompop, &topop);
	    if(frompop==topop)
	      {
		symbol_Theta(lx-40, ly+125, 12, frompop+1);
	      }
	    else
	      {
		symbol_M(lx-40, ly+125, 12, frompop+1, topop+1, world->options->usem);
	      }
	  }
	else
	  {
	    if (i0==world->numpop2 && mu)
	      continue;
	    else
	      {
		s = get_which_species_model(i,world->species_model,world->species_model_size);
		if (i==s->paramindex_mu)
		  {
		    frompop = s->from;
		    topop   = s->to;
		    symbol_D(lx-40, ly+125, 12, frompop+1, topop+1);
		  }
	      }
	  }
	pdf_print_contents_at(lx+90,ly-25, "Time [scaled by mutation rate / generation (DNA: per site, other: per locus)]");
        HPDF_Page_SetFontAndSize (page, helv, 10.0);
        
        // focus on all recorded times or only on 0..mean of largest population
        if(!enlarged)
        {
            // all records
            upperlimit = quantiler(eventbins_all[i], 1.00, 1.00, (long) bins/2,bins);
            if(upperlimit < ((double) (2*world->bayes->histogram[loc].modes[i])))
	      upperlimit =  (double) (2*world->bayes->histogram[loc].modes[i]);
            pdf_histogram_plus(eventbins_all[i], std[i], set50, set95, bins, delta,
                               0., -9999, lx, ly, page_width - 55 - lx, 116, upperlimit, TRUE, world, confidence[i], topop);
        }
        else
        {
            // up to mean of largest population
            loc = (loci > 1) ? loci : 0;
            //            m2mm(i, world->numpop, &frompop, &topop);
            maxtime = (double) world->bayes->histogram[loc].means[0];
            for(pop=1;pop<world->numpop;pop++)
            {
                if(world->bayes->histogram[loc].means[pop]>maxtime)
		  maxtime= (double) world->bayes->histogram[loc].means[pop];
            }
            //xcode j=eventbinnum[i]-1;
            adjustbins = (long) (ceil(maxtime/delta));
            // why could bin be smaller than adjustbins?
            bins = bins < adjustbins ? bins : adjustbins;
            if(bins==0)
                continue;
            lasttimebin = delta * bins;
            // evaluates how much of the y-axis should be shown
            upperlimit = quantiler(eventbins_all[i], 1.0, 0.95, (long) bins/2,bins);
            pdf_histogram_plus(eventbins_all[i], std[i], set50, set95, bins, delta,
                               0., lasttimebin,
                               lx, ly, page_width - 55 - lx, 116, upperlimit ,TRUE, world, confidence[i], topop);
        }
        page_height -= 160;
        if(i < (numparams-1))
        {
            pdf_page_advance_or_not(&page_height, 50);
            ph = (double) HPDF_Page_GetHeight(page) - 55.0 - LINESTRETCH;
            if(page_height >= ph)
                page_height -= 150;
            ly = page_height;
        }
    }
    myfree(binning);
    myfree(set50);
    myfree(set95);
    free_doublevec2d(eventbins_all);
    free_doublevec2d(std);
    free_doublevec2d(confidence);
    myfree(eventbin1max);
}
///
/// Data set was subsampled and used a random sample of size
void pdf_print_random_subset(data_fmt * data, option_fmt *options)
{
    long locus;
    long pop;
    char *name;
    char **header;
    char ***elements;
    double page_width;
    long elements_alloc = data->loci * data->numpop;
    long linenumber;
    long maxnum;
    long count = 0;
    long ind;
    long index;
    long length;
    if(options->randomsubset > 0)
    {
        //print title
      page_width = (double) HPDF_Page_GetWidth(page);
        pdf_print_section_title(&page_width, &page_height, "Subsampled dataset");
        pdf_advance(&page_height);
        pdf_printf(55.,page_height, 'L', "Data set was subsampled randomly per population: %li samples taken",
                   options->randomsubset);
        pdf_advance(&page_height);
        charvec2d(&header, 3,LINESIZE);
        sprintf(header[0],"Locus");
        sprintf(header[1],"Population");
        sprintf(header[2],"Individuals");
        name = (char*) mycalloc(options->nmlength+1,sizeof(char));
        elements = (char ***) mycalloc(elements_alloc, sizeof(char**));
        linenumber=0;
        for (locus=0; locus < data->loci; locus++)
        {
            for (pop=0; pop < data->numpop; pop++)
            {
                if(linenumber>=elements_alloc)
                {
                    elements_alloc += data->loci * data->numpop;
                    elements = (char***) realloc(elements,sizeof(char **) * (size_t) elements_alloc);
                }
                elements[linenumber] = (char **) mycalloc(3,sizeof(char *));
                elements[linenumber][0] = (char *) mycalloc(7, sizeof(char));
                elements[linenumber][1] = (char *) mycalloc(12, sizeof(char));
                elements[linenumber][2] = (char *) mycalloc(LINELENGTH, sizeof(char));
                if(pop==0)
                    sprintf(elements[linenumber][0],"%5li ",locus+1);
                else
                    sprintf(elements[linenumber][0]," ");
                sprintf(elements[linenumber][1], "%-10.10s ", data->popnames[pop]);
                maxnum = options->randomsubset < data->numind[pop][locus] ? options->randomsubset : data->numind[pop][locus];
                count = 0;//18 characters are already consumed on line, see below
                for(ind=0;ind<maxnum;ind++)
                {
                    index = data->shuffled[pop][locus][ind];
                    if(data->indnames[pop][index][locus][0]=='\0')
		      memcpy(name,data->indnames[pop][index][0],sizeof(char)*(size_t) options->nmlength);
                    else
		      memcpy(name,data->indnames[pop][index][locus],sizeof(char)* (size_t) options->nmlength);
                    remove_trailing_blanks(&name);
                    length = (long) strlen(name);
                    if (count+length < LINELENGTH-18)
                    {
                        count += sprintf(elements[linenumber][2] + count ,"%s ",name);
                    }
                    else
                    {
                        sprintf(elements[linenumber][2] + count ,"\n");
                        count = 0;
                        linenumber++;
                        if(linenumber>=elements_alloc)
                        {
                            elements_alloc += data->loci * data->numpop;
                            elements = (char***) realloc(elements,sizeof(char **) * (size_t) elements_alloc);
                        }
                        elements[linenumber] = (char **) mycalloc(3,sizeof(char *));
                        elements[linenumber][0] = (char *) mycalloc(7, sizeof(char));
                        elements[linenumber][1] = (char *) mycalloc(12, sizeof(char));
                        elements[linenumber][2] = (char *) mycalloc(LINELENGTH, sizeof(char));
                        sprintf(elements[linenumber][0]," ");
                        sprintf(elements[linenumber][1]," ");
                        count += sprintf(elements[linenumber][2] + count ,"%s ",name);
                    }
                }
                linenumber++;
            }
        }
        pdf_table2(3, (int) linenumber, header, NULL, elements, NULL, 2, 10.0);
        myfree(name);
        myfree(header);
        for (locus=0; locus < linenumber;locus++)
        {
            myfree(elements[locus][0]);
            myfree(elements[locus][1]);
            myfree(elements[locus][2]);
            myfree(elements[locus]);
        }
        myfree(elements);
    }
}

void pdf_print_spectra(world_fmt *world, data_fmt *data, option_fmt *options, MYREAL ***freq, MYREAL ** total, MYREAL *grandtotal, MYREAL *avghet, MYREAL avghet1, long * maxalleles)
{
    char **header;
    char ***elements;
    long mostalleles =0;
    long pop;
    long locus;
    long a;
    long *maxallelepop;
    double page_width;
    MYREAL allfreq;
    long   sublocus;
    long   sublocistart;
    long   sublociend;
    mutationmodel_fmt *s;
    
    long pop1;
    double homo;
    double f;
    double fx;
    double general_homo;
    double v;

    /* for large datasets this allows to not get overwhelmed with the PDF*/
    if (options->tersepdf)
      return;

    //double numdiv = options->newpops_numpop;
    maxallelepop = (long *) mycalloc(data->numpop, sizeof(long));
    charvec2d(&header, data->numpop+2,LINESIZE);
    //print title
    page_width = (double) HPDF_Page_GetWidth(page);
    pdf_print_section_title(&page_width, &page_height, "Allele frequency spectra");
    // loop over loci
    pdf_advance(&page_height);
    // header
    sprintf(header[0],"Allele");
    for(pop1=0; pop1 < data->numpop; pop1++)
    {
        pop = options->newpops[pop1]-1;
        sprintf(header[pop1+1],"Pop%-3li",pop+1);
    }
    sprintf(header[pop1+1],"All");
    for(locus=0; locus < data->allsubloci; locus++)
    {
        s= &world->mutationmodels[locus];
        if(mostalleles < s->numstates)
            mostalleles = s->numstates;
        //if(mostalleles < s->maxalleles)
	//  mostalleles = s->maxalleles;
    }
    // max of alleles + Total + H_exp; this is done for each locus
    elements = (char ***) mycalloc(mostalleles+3, sizeof(char**));
    for(a=0; a < mostalleles+3; a++)
    {
        charvec2d(&elements[a],data->numpop+2, LINESIZE);
    }
    
    for(locus=0; locus < data->loci; locus++)
    {
      memset(maxallelepop,0,sizeof(long)*(size_t) data->numpop);
        pdf_printf(55.,page_height, 'L', "Locus %i",locus+1);
        pdf_advance(&page_height);
        sublocistart = world->sublocistarts[locus];
        sublociend = world->sublocistarts[locus+1];
        general_homo=0.0;
        for(sublocus=sublocistart;sublocus<sublociend;sublocus++)
        {
            s= &world->mutationmodels[sublocus];
            //for(a=0; a < s->maxalleles; a++)
	    for(a=0; a < s->numstates; a++)
            {
                sprintf(elements[a][0],"%s ",data->allele[sublocus][a]);
                allfreq = 0.0;
                for(pop1=0; pop1 < data->numpop; pop1++)
                {
                    pop = options->newpops[pop1]-1;
                    // freq was constructed using newpop!!
                    if (freq[pop][locus][a]>0.0)
                    {
                        maxallelepop[pop1] += 1;
                        sprintf(elements[a][pop1+1],"%1.3f",freq[pop][sublocus][a]/total[pop][sublocus]);
                        allfreq += freq[pop][sublocus][a];
                    }
                    else
                    {
                        sprintf(elements[a][pop1+1],"  -  ");
                    }
                }
                sprintf(elements[a][pop1+1],"%1.3f",allfreq/grandtotal[sublocus]);
                fx  =  (double) (allfreq/grandtotal[sublocus]);
                general_homo += fx * fx;
            }
            sprintf(elements[a][0],"Alleles");
            for(pop1=0; pop1 < data->numpop; pop1++)
            {
                sprintf(elements[a][pop1+1],"%li",maxallelepop[pop1]);
            }
            //sprintf(elements[a][pop1+1],"%li",s->maxalleles);
	    sprintf(elements[a][pop1+1],"%li",s->numstates);
	    sprintf(elements[a+1][0],"Samplesize");
	    for(pop1=0; pop1 < data->numpop; pop1++)
	      {
		sprintf(elements[a+1][pop1+1],"%li", (long) total[pop1][locus]);
	      }
	    sprintf(elements[a+1][pop1+1],"%li",(long) grandtotal[locus]);
            sprintf(elements[a+2][0],"_H_exp");
            for (pop1 = 0; pop1 < data->numpop; pop1++)
            {
                pop = options->newpops[pop1]-1;
                homo = 0.0;
		long aa;
                for(aa=0;aa<maxalleles[locus];aa++)
                {
		  f = (double) (freq[pop][locus][aa]/total[pop][sublocus]);
                    homo += f*f;
                }
                v = 1.0 - homo;
                //avghet1 += v;
                //avghet[pop1] += v;
                sprintf(elements[a+2][pop1+1],"%5.3f",v);
            }
            sprintf(elements[a+2][pop1+1]," %5.3f",1.0-general_homo);
            //avghet1 += 1.0 - general_homo;
            //pdf_table2((int) (data->numpop+2),(int) (s->maxalleles+3), header, NULL, elements, NULL, 2, 10.0);
	    pdf_table2((int) (data->numpop+2),(int) (s->numstates+3), header, NULL, elements, NULL, 2, 10.0);
            pdf_advance(&page_height);
        }
    }
    pdf_printf(55.,page_height, 'L', "Average expected heterozygosity");
    pdf_advance(&page_height);
    for (pop1 = 0; pop1 < data->numpop; pop1++)
    {
      //pop = options->newpops[pop1]-1;
        sprintf(elements[0][pop1+1], "%5.3f ",avghet[pop1] / data->loci);
        //avghetall += avghet[pop1] / data->loci;
    }
    sprintf(elements[0][pop1+1],"%5.3f",avghet1/data->loci);
    strcpy(header[0]," ");
    strcpy(elements[0][0],"_H_exp");
    pdf_table2((int)(data->numpop+2),1, header, NULL, elements, NULL, 2, 10.0);
    free_charvec2d(header);
    for(a=0; a < mostalleles+3; a++)
        free_charvec2d(elements[a]);
    myfree(elements);
    myfree(maxallelepop);
}

///
/// print average temperatures
void pdf_print_averageheat(world_fmt **universe, option_fmt *options)
{
    long a;
    char **header;
    char ***elements;
    //  double position[10] = {55., 100., 150., 200., 250., 300., 350., 400., 450., 500.};
    double page_width;
    charvec2d(&header, 2,LINESIZE);
    //print title
    page_width = (double) HPDF_Page_GetWidth(page);
    pdf_print_section_title(&page_width, &page_height, "Average temperatures during the run");
    pdf_advance(&page_height);
    // header
    sprintf(header[0],"Chain");
    sprintf(header[1],"Temperatures");
    elements = (char ***) mycalloc(options->heated_chains, sizeof(char**));
    for(a=0; a < options->heated_chains; a++)
        charvec2d(&elements[a],2, LINESIZE);
    for(a=0; a < options->heated_chains; a++)
    {
        sprintf(elements[a][0],"%5li ",a+1);
	if (options->adaptiveheat == STANDARD)
	  sprintf(elements[a][1],"%10.5f ",universe[a]->averageheat);
	else
	  sprintf(elements[a][1],"%10.5f ",universe[a]->heat);
    }
    pdf_table2( 2, (int) (options->heated_chains), header, NULL, elements, NULL, 2, 10.0);
    pdf_advance(&page_height);
    if (options->adaptiveheat == STANDARD)
      {
	pdf_advance(&page_height);
	pdf_printf_next(55., &page_height,"Adaptive heating often fails, if the average temperatures are very close together\n");
	pdf_printf_next(55., &page_height,"try to rerun using static heating! If you want to compare models using marginal\n");
	pdf_printf_next(55., &page_height,"likelihoods then you MUST use static heating\n");
	pdf_advance(&page_height);
      }
    free_charvec2d(header);
    for(a=0; a < options->heated_chains; a++)
        free_charvec2d(elements[a]);
    myfree(elements);
}

///
/// table frequency of events through time for each locus and all loci
void pdf_print_eventtime_table(world_fmt *world)
{
    double   page_width;
    double   w;
    long    i;
    long    j;
    long    j0;
    long    frompop;
    long    topop;
    long   *eventbinnum;
    long    locus;
    duo   **eventbins;
    double   age;
    char    title[LINESIZE];
    //double   left_margin = 55;
    long     end = world->loci > 1 ? world->loci + 1 : 1;
    long     start = 0;
    sprintf(title,"Distribution of events trough time");
    pdf_new_page("");
    HPDF_Page_SetFontAndSize (page, helvob, 16.0);
    w = (double) HPDF_Page_TextWidth(page, title);
    
    /* Start to print text. */
    page_height = (double) HPDF_Page_GetHeight(page);
    page_width = (double) HPDF_Page_GetWidth(page);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    HPDF_Page_SetRGBStroke(page, 0,0,0);
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);

    /* for large datasets this allows to not get overwhelmed with the PDF*/
    if (world->options->tersepdf)
      {       
	start = end - 1;
      }
    for(locus=start;locus < end; locus++)
    {
        if(locus<world->loci)
        {
            if(world->data->skiploci[locus])
                continue;
        }
        eventbins = world->mighistloci[locus].migeventbins;
        eventbinnum = world->mighistloci[locus].migeventbinnum;
        // population sizes
        for(j0=0; j0 < world->numpop2 + world->bayes->mu + 2 * world->species_model_size; j0++)
        {
            if(!shortcut(j0,world,&j))
            {
                pdf_advance(&page_height);
                if(locus == world->loci)
                    pdf_printf(left_margin, page_height, 'L' , "All loci");
                else
                    pdf_printf(left_margin, page_height, 'L' , "Locus %li:",locus);
                if(j0 < world->numpop)
                    symbol_Theta(left_margin+80, page_height, 12, j0+1);
                else
                {
		  if(j < world->numpop2)
		    {
		      m2mm (j0, world->numpop, &frompop, &topop);
		      symbol_M(left_margin+80, page_height, 12, frompop+1, topop+1, world->options->usem);
		    }
		  else
		    {
		      if(world->bayes->mu && j==world->numpop2)
			continue; //rate		      
		      species_fmt * s = get_which_species_model(j,world->species_model,world->species_model_size);
		      if ( j == s->paramindex_mu)
			symbol_D(left_margin+80, page_height, 12, s->from+1, s->to+1);
		      else
			continue;
		    }
		}
		pdf_advance(&page_height);
		pdf_printf(left_margin, page_height, 'L', "Time");
                pdf_printf(left_margin + 200, page_height, 'L', "Frequency of visit");
                pdf_printf(left_margin + 400, page_height, 'L', "MRCA frequency");
                pdf_advance_half(&page_height);
                pdf_draw_line(50, page_height, page_width-50, page_height);
                pdf_advance(&page_height);
                age = (double) (-world->mighistloci[locus].eventbinsize/2.);
                for(i = 0; i < eventbinnum[j]; i++)
                {
                    age += world->mighistloci[locus].eventbinsize;
                    pdf_printf(left_margin, page_height, 'L', "%10.10f", age);
                    pdf_printf(left_margin+200, page_height, 'L', "%8.5f", (MYREAL) eventbins[j][i][0]);
                    pdf_printf(left_margin+400, page_height, 'L', "%8.5f", (MYREAL) eventbins[j][i][1]);
                    pdf_advance(&page_height);
                }
                pdf_draw_line(50, page_height, page_width-50, page_height);
                pdf_advance(&page_height);
            }
        }
    }
}
///
/// average and median time for a events
void pdf_print_time_table(world_fmt *world,
                          double *meantime, double *mediantime, double *stdtime, double *freq,
                          boolean mrca)
{
    double   page_width;
    double   w;
    long    pop;
    long pop0;
    long    lp;
    long    frompop;
    long    topop;
    long    locus;
    long    npall = world->numpop2+world->species_model_size * 2 + world->bayes->mu;
    char    title[LINESIZE];
    //double   left_margin = 55;
    long     end = world->loci > 1 ? world->loci + 1: 1;
    long     start = 0;

    if(mrca)
        sprintf(title,"Time and probability of location of most recent common ancestor");
    else
        sprintf(title,"Summary statistics of events through time");
    pdf_new_page("");
    HPDF_Page_SetFontAndSize (page, helvob, 16.0);
    w = (double) HPDF_Page_TextWidth(page, title);
    
    /* Start to print text. */
    page_height = (double) HPDF_Page_GetHeight(page);
    page_width = (double) HPDF_Page_GetWidth(page);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    HPDF_Page_SetRGBStroke(page, 0,0,0 );
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
    if (world->options->tersepdf)
      {
	start = end - 1;
      }
    for(locus=start;locus < end; locus++)
    {
        if(world->data->skiploci[locus])
            continue;
        
        //xcode eventbins = world->mighistloci[locus].migeventbins;
        //xcode eventbinnum = world->mighistloci[locus].migeventbinnum;
        
        pdf_advance(&page_height);
        if(locus == world->loci)
            pdf_printf(left_margin, page_height, 'L' , "All loci");
        else
            pdf_printf(left_margin, page_height, 'L' , "Locus %li",locus+1);
        
        pdf_advance(&page_height);
        pdf_printf(left_margin, page_height, 'L', "Population");
        pdf_printf(left_margin + 100, page_height, 'L', "Time");
        pdf_printf(left_margin + 400, page_height, 'L', "Frequency");
        pdf_advance_half(&page_height);
        pdf_draw_line(left_margin + 100, page_height, left_margin+390, page_height);
        pdf_advance(&page_height);
        if(!mrca)
        {
            pdf_printf(left_margin, page_height, 'L', "From");
            pdf_printf(left_margin + 50, page_height, 'L', "To");
        }
        pdf_printf(left_margin + 100, page_height, 'L', "Average");
        pdf_printf(left_margin + 200, page_height, 'L', "Median");
        pdf_printf(left_margin + 300, page_height, 'L', "Std");
        pdf_advance_half(&page_height);
        pdf_draw_line(50, page_height, page_width-50, page_height);
        pdf_advance(&page_height);
        for (pop0 = (world->options->mighist_all ? 0 : world->numpop);
             pop0 <  (mrca ? world->numpop : npall); pop0++)
	  {
	    if(shortcut(pop0, world,&pop))
	      {
		continue;
	      }
	    if(pop<world->numpop2)
	      m2mm(pop, world->numpop, &frompop, &topop);
	    else
	      {
		if (pop==world->numpop2 && world->bayes->mu)
		  {
		    warning("Table of events for rates is not implemented");
		    continue;
		  }
		else
		  {
		    species_fmt * s = get_which_species_model(pop,world->species_model, world->species_model_size);
		    if (pop == s->paramindex_mu)
		      {
			frompop = s->from;
			topop = s->to;
		      }
		    else
		      continue;
		  }
	      }

            lp = locus * world->numpop2 + pop;
            if(freq[lp]>0.0)
            {
                pdf_printf(left_margin,       page_height, 'L', "%li",frompop+1);
                if(!mrca)
                    pdf_printf(left_margin + 50,  page_height, 'L', "%li", topop+1);
                pdf_printf(left_margin + 100, page_height, 'L', "%f", meantime[lp]);
                pdf_printf(left_margin + 200, page_height, 'L', "%f", mediantime[lp]);
                pdf_printf(left_margin + 300, page_height, 'L', "%f", stdtime[lp]);
                pdf_printf(left_margin + 400, page_height, 'L', "%f", freq[lp]);
                pdf_advance(&page_height);
            }
        }
        pdf_advance(&page_height);
    }
}


///
/// plot event histograms
void
pdf_histogram_legend(void)
{
    double   page_width;
    double   lx;
    double   w;
    char    title[LINESIZE];
    
    // add a new page for the legend
    pdf_new_page("");
    sprintf(title,"%s","Legend for Skyline and Event plots");
    HPDF_Page_SetFontAndSize (page, helvob, 16.0);
    w = (double) HPDF_Page_TextWidth(page, title);
    
    /* Start to print text. */
    page_height = (double) HPDF_Page_GetHeight(page);
    page_width = (double) HPDF_Page_GetWidth(page);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    HPDF_Page_SetRGBStroke(page, 0.f, 0.f, 0.f);
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    HPDF_Page_SetFontAndSize (page, helv, 10.0);
    // first histogram position
    pdf_advance(&page_height);
    lx = 60;
    pdf_printf_next(lx-5, &page_height, "Skyline plots:");
    pdf_printf_next(lx, &page_height, "Skyline plots visualize the changes of population sizes and migration rates through time");
    pdf_printf_next(lx, &page_height, "(today is on the left side and time is measured into the past. The time scale is in units of");
    pdf_printf_next(lx, &page_height, "expected mutations per generation. To calculate the absolute time scale you must supply an");
    pdf_printf_next(lx, &page_height, "mutation rate per year and the duration of a  generation in years in the data option.");
    pdf_printf_next(lx, &page_height, "You can calculate the absolute time by multiplying the scale by generation time times ");
    pdf_printf_next(lx, &page_height, "mutation rate per year (per site for DNA; per locus for all other datatypes).");
    pdf_advance(&page_height);
    pdf_printf_next(lx, &page_height, "With estimated mutation rate only the combined rate modifier is plotted.");
    pdf_printf_next(lx, &page_height, "[this will change to  mutation rate plot].");
    pdf_advance(&page_height);
    pdf_printf_next(lx, &page_height, "The gray bars cover 1.96 * approximate standard error (std in file skylinefile/number");
    pdf_printf_next(lx, &page_height, "of observations in the bin) up and down from the expected value.");
    pdf_printf_next(lx, &page_height, "The bar with different shades of gray on top of each plot indicates the number of values that were used to");
    pdf_printf_next(lx, &page_height, "to calculate the expected value, white means there were very few and black means");
    pdf_printf_next(lx, &page_height, "that there were many thousands of samples per bin.");
    pdf_advance(&page_height);
    pdf_printf_next(lx, &page_height, "On some plots one can see red squares below the grayscale bar, these suggest that either the");
    pdf_printf_next(lx, &page_height, "upper quantile and/or the main value was higher than the visible part of the  axis.");
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    
    pdf_printf_next(lx-5, &page_height, "Event histograms:");
    pdf_printf_next(lx, &page_height, "All accepted events (migration events, coalescent events) are recorded and their frequency");
    pdf_printf_next(lx, &page_height, "are shown as histograms over time with recent time on the left side. The frequency plots of");
    pdf_printf_next(lx, &page_height, "populations with constant size and constant immigration rates show histograms that are similar");
    pdf_printf_next(lx, &page_height, "to exponential distribution, if the populations come from a divergence model without migration");
    pdf_printf_next(lx, &page_height, "then the frequency of migration events can show a peak in the past.");
}

void pdf_print_haplotype_title(void)
{
    //  page_height = HPDF_Page_GetHeight(page);
    double page_width = (double) HPDF_Page_GetWidth(page);
    pdf_print_section_title(&page_width, &page_height, "Haplotype states and probabilities");
}

void pdf_print_endline(void)
{
    pdf_advance(&page_height);
}
//
// print haplotype states
void pdf_print_haplotypes2(char *buffer, long buflen)
{
    
  //double left_margin  = 55.;
    double right_margin = (double) HPDF_Page_GetWidth(page) - 55.0;
    long b=0;
    double x=left_margin;
    double namewidth = (double) HPDF_Page_TextWidth(page, "0123456789  ");
    for(b=0;b<buflen;b++)
    {
        char bufchar = buffer[b];
        switch(bufchar)
        {
            case '\t':
                x = left_margin + namewidth;
                break;
            case '\n':
                if(buffer[b+1]=='\t')
                    x = left_margin + namewidth;
                else
                    x = left_margin;
                pdf_advance(&page_height);
                break;
            default:
                pdf_putc(&x, &page_height, left_margin, right_margin, bufchar);
                break;
        }
    }
}

/*
 void pdf_print_haplotypes(world_fmt *world)
 {
 long locus;
 long ind;
 
 //  for(locus=0;locus<world->loci;locus++)
 //  {
 //    for(ind=0; ind<world->data->haplotyping_numind[locus];ind++)
 //	{
 pdf_print_haplotypes2(world->data->haplotyping_buffer[locus][ind],world->data->haplotyping_buflen[locus][ind]);
 //	}
 // }
 }
 */


#ifdef BEAGLE
void print_beagle_resources_pretty(BeagleInstanceDetails instDetails)
{
    pdf_advance(&page_height);
    pdf_print_contents_at(left_margin, page_height,"Likelihood calculator using HMSBEAGLE:");
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height, 'L', "\tRsrc Name : %s", instDetails.resourceName);
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height, 'L', "\tImpl Name : %s\n", instDetails.implName);
    pdf_advance(&page_height);
}
#endif

///
///
void findminmax(double *vals, const long n, double *min, double *max)
{
    long i;
#ifdef WINDOWS
    *min = 1000000;
    *max = -1000000;
#else
    *min = (double) HUGE;
    *max = (double) -HUGE;
#endif
    for(i=0; i < n; i++)
    {
        double v = (double) vals[i];
        if(*min > v)
            *min = v;
        if(*max < v)
            *max = v;
    }
    //  printf("min=%f, max=%f\n",*min,*max);
}

///
/// Line plot assumes that x and y are a series of coordinates, connection with a colored line with given thickness, dots are marked in color
void pdf_linedotplot(long n, double *x, double *y, double dotcolor[], double linecolor[], MYREAL linethickness, double width, double height, boolean has_dots, boolean has_line)
{
    double minx;
    double miny;
    double maxx;
    double maxy;
    //  double page_width = HPDF_Page_GetWidth(page);
    double lx = 100;
    double ly = page_height - height;
    long i;
    float red = (float) linecolor[0];
    float green = (float) linecolor[1];
    float blue = (float) linecolor[2];
    double xrange;
    double yrange;
    page_height = ly;
    pdf_advance(&page_height);
    findminmax(x, n, &minx, &maxx);
    findminmax(y, n, &miny, &maxy);
    pdf_create_axes(VERTICAL, miny, maxy, 5, lx, ly, height);
    pdf_create_axes(HORIZONTAL, minx , maxx, 9, lx, ly, width);
    if(has_line)
    {
      HPDF_Page_SetLineWidth(page, (float) linethickness);
        HPDF_Page_SetRGBStroke(page, red, green, blue);
        xrange = maxx - minx;
        yrange = maxy - miny;
        width /= xrange;
        height/= yrange;
        for(i=0; i < n-1; i++)
        {
            pdf_draw_line(lx+(x[i]-minx)*width, ly+(y[i]-miny)*height,lx+(x[i+1]-minx)*width,ly+(y[i+1]-miny)*height);
        }
        HPDF_Page_SetRGBStroke(page, 0.f,0.f,0.f);//reset to black
    }
    if(has_dots)
    {
        for(i=0; i < n; i++)
        {
            pdf_print_dot(lx+(x[i]-minx)*width, ly+(y[i]-miny)*height, 3, SQUARE, dotcolor);
        }
    }
}


/*
 ==============================================================================
 Likelihood ratio tests
 ==============================================================================
 Over all loci
 Legend for the LRT tables
 -------------------------------------------------------------------------------
 Null-Hypothesis: your test model         | Log(likelihood) of test model
 =same=                                   | Log(likelihood) of full model
 full model (the model under which the    | Likelihood ratio test value
 genealogies were sampled)                | Degrees of freedom of test
 [Theta values are on the diagonal of the | Probability*
 Migration matrix, migration rates are    | Probability**
 specified as M]                          | Akaike's Information Criterion***
 | Number of parameters used
 -------------------------------------------------------------------------------
 *) Probability under the assumption that parameters have range -Inf to Inf
 **) Probability under the assumption that parameters have range 0 to Inf
 ***) AIC: the smaller the value the better the model
 [the full model has AIC=-2.396909, num(param)=4]
 */
//#define OFFSET 300
//#define OFFSET2 280
//#define ADDED 70

// citations
void pdf_print_citation(char type[],world_fmt *world)
{
    char mytype = type[0];
    pdf_printf(left_margin, page_height,'L',"%s","Citation suggestions:"); 
    pdf_advance(&page_height);
    switch(mytype)
    {
        case 'M': // marginal likelihood
            pdf_printf(left_margin, page_height,'L',"Beerli P. and M. Palczewski, 2010. Unified framework to evaluate panmixia and migration direction among");
            pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"    multiple sampling locations, Genetics, 185: 313-326.");
            pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"Palczewski M. and P. Beerli, 2014. Population model comparison using multi-locus datasets.");
	    pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"    In M.-H. Chen, L. Kuo, and P. O. Lewis, editors, Bayesian Phylogenetics: Methods, ");
	    pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"    Algorithms, and Applications, pages 187-200. CRC Press, 2014.");
            //pdf_advance(&page_height);
            //pdf_printf(left_margin, page_height,'L',"Xie W., P. O. Lewis, Y. Fan, L. Kuo, and M.-H. Chen. 2011. Improving marginal likelihood");
	    //pdf_advance(&page_height);
            //pdf_printf(left_margin, page_height,'L',"    estimation for Bayesian phylogenetic model selection. Systematic Biology, 60(2):150–160, 2011.");
	    pdf_advance(&page_height);
            break;
        case 'L': // maximum likelihood
            pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"Beerli P., 1998. Estimation of migration rates and population sizes in geographically structured populations.");
            pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"    In Advances in Molecular Ecology, G. R. Carvalho, ed., vol. 306 of NATO sciences series, Series A: Life sciences,");
            pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"    ISO Press, Amsterdam, pp. 39-53.");
            pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"Beerli P. and J. Felsenstein, 1999. Maximum-likelihood estimation of migration rates and effective population");
            pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"    numbers in two populations using a coalescent approach, Genetics, 152:763-773.");
            pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"Beerli P. and J. Felsenstein, 2001. Maximum likelihood estimation of a migration matrix and effective ");
            pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"    population sizes in n subpopulations by using a coalescent approach, Proceedings of the National Academy ");
            pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"    of Sciences of the United States of America, 98: p. 4563-4568.");
            pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"Beerli P., 2007. Estimation of the population scaled mutation rate from microsatellite data,");
            pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"    Genetics, 177:1967-1968.");
            pdf_advance(&page_height);
            if (world->options->datatype=='m' || world->options->datatype=='s')
            {
                pdf_printf(left_margin, page_height,'L',"Beerli P., 2007. Estimation of the population scaled mutation rate from microsatellite data,");
                pdf_advance(&page_height);
                pdf_printf(left_margin, page_height,'L',"    Genetics, 177:1967-1968.");
                pdf_advance(&page_height);
            }
        case 'B': // Bayesian inference
            pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"Beerli P., 2006. Comparison of Bayesian and maximum-likelihood inference of population genetic parameters.");
            pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"    Bioinformatics 22:341-345");
            pdf_advance(&page_height);

	    pdf_printf(left_margin, page_height,'L',"Beerli, P., H. Ashki, S. Mashayekhi, and M. Palczewski, 2022. Population divergence time estimation");
	    pdf_advance(&page_height);
	    pdf_printf(left_margin, page_height,'L',"   using individual lineage label switching. G3 Genes—Genomes—Genetics, 12(4), 02 2022.");
            pdf_advance(&page_height);
	    
            if (world->options->datatype=='m' || world->options->datatype=='s')
            {
                pdf_printf(left_margin, page_height,'L',"Beerli P., 2007. Estimation of the population scaled mutation rate from microsatellite data,");
                pdf_advance(&page_height);
                pdf_printf(left_margin, page_height,'L',"    Genetics, 177:1967-1968.");
                pdf_advance(&page_height);
            }
            pdf_printf(left_margin, page_height,'L',"Beerli P., 2009. How to use MIGRATE or why are Markov chain Monte Carlo programs difficult to use?");
            pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"    In Population Genetics for Animal Conservation, G. Bertorelle, M. W. Bruford, H. C. Hauffe, A. Rizzoli,");
            pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"    and C. Vernesi, eds., vol. 17 of Conservation Biology, Cambridge University Press, Cambridge UK, pp. 42-79.");
            pdf_advance(&page_height);

	    pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"Beerli, P., S. Mashayekhi, M. Sadeghi, M. Khodaei, and K. Shaw, 2019. Population genetic inference with migrate.");
	    pdf_advance(&page_height);
            pdf_printf(left_margin, page_height,'L',"    Current Protocols in Bioinformatics, 68(1):e87.");
	    pdf_advance(&page_height);
    }
}

void pdf_report_unassigned(world_fmt *world)
{
    long i;
    long idi;
    long locus;
    long pop;
    MYREAL sum = 0.0;
    MYREAL lsum;
    MYREAL *total;
    MYREAL totalsum;
    char *key;
    char **header;
    char **header2;
    long a;
    double w;
    char ***elements;
    char ***elements_sum;
    char *title;
    double page_width;
    long siz;
    long siz_sum;
    long ii;
    long z;
    if (world->has_unassigned)
    {
        siz_sum = world->unassignednum;
        siz = world->unassignednum*(world->loci+1);
        total = (MYREAL *) mycalloc(world->numpop, sizeof(MYREAL));
        charvec2d(&header, 2* (2 + world->numpop),STRSIZE);
        header2 = header + 2 + world->numpop;
        elements = (char ***) mycalloc(siz, sizeof(char**));
        for(a=0; a<siz;a++)
            charvec2d(&elements[a],  2 + world->numpop, STRSIZE );
        elements_sum = (char ***) mycalloc(siz_sum, sizeof(char**));
        for(a=0; a<siz_sum;a++)
            charvec2d(&elements_sum[a],  2 + world->numpop, STRSIZE );
        sprintf(header[0],"Individual");
        sprintf(header[2],"Population");
        /*      for (pop=0;pop<world->numpop;pop++)
         {
         sprintf(header[pop+2],"     ");
         }*/
        for (pop=0;pop<world->numpop;pop++)
        {
            sprintf(header2[2+pop],"%5li",pop+1);
        }
        ii=0;
        z=0;
        for (i=1;i<world->unassignednum;i++)
        {
            key = world->unassigned[i]->key;
            memset(total,0,sizeof(MYREAL)*(size_t) world->numpop);
            for (locus=0;locus<world->loci;locus++)
            {
                sprintf(elements[ii*(world->loci+1)+locus][0], "%-10.10s",key);
                sum = 0.0;
                for (pop=0;pop<world->numpop;pop++)
                {
                    sum += world->unassigned[i]->probloc[INDIX(world->numpop,locus,pop)];
                }
                sprintf(elements[ii*(world->loci+1)+locus][1], "%8li",locus+1);
                for (pop=0;pop<world->numpop;pop++)
                {
                    idi = INDIX(world->numpop,locus,pop);
                    sprintf(elements[ii*(world->loci+1)+locus][pop+2], "%5.3f", world->unassigned[i]->probloc[idi]/sum);
                    total[pop] += log(world->unassigned[i]->probloc[idi]/sum);	          
                }
            }
            totalsum = 0.0;
            for (pop=0;pop<world->numpop;pop++)
            {
                totalsum += exp(total[pop]);
            }
            lsum = log(totalsum);
            sprintf(elements[ii*(world->loci+1)+locus][0], "%-10.10s",key);
            sprintf(elements[ii*(world->loci+1)+locus][1], "     All  ");
            sprintf(elements_sum[z][0], "%-10.10s",key);
            sprintf(elements_sum[z][1], "          ");
            for (pop=0;pop<world->numpop;pop++)
            {
                sprintf(elements[ii*(world->loci+1)+locus][pop+2], "%5.3f", exp(total[pop]-lsum));
                sprintf(elements_sum[z][pop+2], "%5.3f", exp(total[pop]-lsum));
            }
            z++;
            ii++;
        }
        // short table
        title = (char *) mycalloc(LINESIZE, sizeof(char));
        sprintf(title, "Summary Assignment of Individuals to Populations");
        pdf_new_page("");
        HPDF_Page_SetFontAndSize (page, helvob, 16.0);
        w = (double) HPDF_Page_TextWidth(page, title);      
        /* Start to print text. */ 
        page_height = (double) HPDF_Page_GetHeight(page);
        page_width = (double) HPDF_Page_GetWidth(page);
        pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
        HPDF_Page_SetRGBStroke(page, 0,0,0);
        page_height -= 126;
        pdf_draw_line(50, page_height, page_width-50, page_height);
        HPDF_Page_SetFontAndSize (page, helv, 10.0);
        pdf_advance(&page_height);
        pdf_advance(&page_height);
        pdf_table2((int)(world->numpop+2), (int) siz_sum, header, header2, elements_sum, NULL, 2, 10.0);      
        //long table
        if (world->loci>1)
        {
            sprintf(title, "Detailed Assignment of Individuals to Populations");
            pdf_new_page("");
            HPDF_Page_SetFontAndSize (page, helvob, 16.0);
            w = (double) HPDF_Page_TextWidth(page, title);      
            page_height = (double) HPDF_Page_GetHeight(page);
            page_width = (double) HPDF_Page_GetWidth(page);
            pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
            HPDF_Page_SetRGBStroke(page, 0,0,0);
            page_height -= 126;
            pdf_draw_line(50, page_height, page_width-50, page_height);
            HPDF_Page_SetFontAndSize (page, helv, 10.0);
            pdf_advance(&page_height);
            pdf_advance(&page_height);
            sprintf(header[1],"Locus     ");
            pdf_table2((int) ( world->numpop+2), (int) siz, header, header2, elements, NULL, 2, 10.0);      
        }
        myfree(total);
        free_charvec2d(header);
        for(a=0; a < siz; a++)
            free_charvec2d(elements[a]);
        myfree(elements);
        for(a=0; a < siz_sum; a++)
            free_charvec2d(elements_sum[a]);
        myfree(elements_sum);
        myfree(title);
    }
}


//#endif /*end of PRETTY*/
