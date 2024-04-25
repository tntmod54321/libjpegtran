/*
 * jpegtran.c
 *
 * This file was part of the Independent JPEG Group's software:
 * Copyright (C) 1995-2019, Thomas G. Lane, Guido Vollbeding.
 * libjpeg-turbo Modifications:
 * Copyright (C) 2010, 2014, 2017, 2019-2022, 2024, D. R. Commander.
 * For conditions of distribution and use, see the accompanying README.ijg
 * file.
 *
 * This file contains a command-line user interface for JPEG transcoding.
 * It is very similar to cjpeg.c, and partly to djpeg.c, but provides
 * lossless transcoding between different JPEG file formats.  It also
 * provides some lossless and sort-of-lossless transformations of JPEG data.
 */

#ifdef _MSC_VER
#define _CRT_SECURE_NO_DEPRECATE
#endif

#include "cdjpeg.h"             /* Common decls for cjpeg/djpeg applications */
#include "transupp.h"           /* Support routines for jpegtran */
#include "jversion.h"           /* for version message */
#include "jconfigint.h"


/*
 * Argument-parsing code.
 * The switch parser is designed to be useful with DOS-style command line
 * syntax, ie, intermixed switches and file names, where only the switches
 * to the left of a given file name affect processing of that file.
 * The main program in this file doesn't actually use this capability...
 */


static const char *progname;    /* program name for error messages */
static char *icc_filename;      /* for -icc switch */
static JDIMENSION max_scans;    /* for -maxscans switch */
static char *outfilename;       /* for -outfile switch */
static char *dropfilename;      /* for -drop switch */
static boolean report;          /* for -report switch */
static boolean strict;          /* for -strict switch */
static JCOPY_OPTION copyoption; /* -copy switch */
static jpeg_transform_info transformoption; /* image transformation options */


LOCAL(void)
usage(void)
/* complain about bad command line */
{
  fprintf(stderr, "usage: %s [switches] ", progname);
#ifdef TWO_FILE_COMMANDLINE
  fprintf(stderr, "inputfile outputfile\n");
#else
  fprintf(stderr, "[inputfile]\n");
#endif

  fprintf(stderr, "Switches (names may be abbreviated):\n");
  fprintf(stderr, "  -copy none     Copy no extra markers from source file\n");
  fprintf(stderr, "  -copy comments Copy only comment markers (default)\n");
  fprintf(stderr, "  -copy icc      Copy only ICC profile markers\n");
  fprintf(stderr, "  -copy all      Copy all extra markers\n");
#ifdef ENTROPY_OPT_SUPPORTED
  fprintf(stderr, "  -optimize      Optimize Huffman table (smaller file, but slow compression)\n");
#endif
#ifdef C_PROGRESSIVE_SUPPORTED
  fprintf(stderr, "  -progressive   Create progressive JPEG file\n");
#endif
  fprintf(stderr, "Switches for modifying the image:\n");
#if TRANSFORMS_SUPPORTED
  fprintf(stderr, "  -crop WxH+X+Y  Crop to a rectangular region\n");
  fprintf(stderr, "  -drop +X+Y filename          Drop (insert) another image\n");
  fprintf(stderr, "  -flip [horizontal|vertical]  Mirror image (left-right or top-bottom)\n");
  fprintf(stderr, "  -grayscale     Reduce to grayscale (omit color data)\n");
  fprintf(stderr, "  -perfect       Fail if there is non-transformable edge blocks\n");
  fprintf(stderr, "  -rotate [90|180|270]         Rotate image (degrees clockwise)\n");
#endif
#if TRANSFORMS_SUPPORTED
  fprintf(stderr, "  -transpose     Transpose image\n");
  fprintf(stderr, "  -transverse    Transverse transpose image\n");
  fprintf(stderr, "  -trim          Drop non-transformable edge blocks\n");
  fprintf(stderr, "                 with -drop: Requantize drop file to match source file\n");
  fprintf(stderr, "  -wipe WxH+X+Y  Wipe (gray out) a rectangular region\n");
#endif
  fprintf(stderr, "Switches for advanced users:\n");
#ifdef C_ARITH_CODING_SUPPORTED
  fprintf(stderr, "  -arithmetic    Use arithmetic coding\n");
#endif
  fprintf(stderr, "  -icc FILE      Embed ICC profile contained in FILE\n");
  fprintf(stderr, "  -restart N     Set restart interval in rows, or in blocks with B\n");
  fprintf(stderr, "  -maxmemory N   Maximum memory to use (in kbytes)\n");
  fprintf(stderr, "  -maxscans N    Maximum number of scans to allow in input file\n");
  fprintf(stderr, "  -outfile name  Specify name for output file\n");
  fprintf(stderr, "  -report        Report transformation progress\n");
  fprintf(stderr, "  -strict        Treat all warnings as fatal\n");
  fprintf(stderr, "  -verbose  or  -debug   Emit debug output\n");
  fprintf(stderr, "  -version       Print version information and exit\n");
  fprintf(stderr, "Switches for wizards:\n");
#ifdef C_MULTISCAN_FILES_SUPPORTED
  fprintf(stderr, "  -scans FILE    Create multi-scan JPEG per script FILE\n");
#endif
  exit(EXIT_FAILURE);
}


LOCAL(void)
select_transform(JXFORM_CODE transform)
/* Silly little routine to detect multiple transform options,
 * which we can't handle.
 */
{
#if TRANSFORMS_SUPPORTED
  if (transformoption.transform == JXFORM_NONE ||
      transformoption.transform == transform) {
    transformoption.transform = transform;
  } else {
    fprintf(stderr, "%s: can only do one image transformation at a time\n",
            progname);
    usage();
  }
#else
  fprintf(stderr, "%s: sorry, image transformation was not compiled\n",
          progname);
  exit(EXIT_FAILURE);
#endif
}


LOCAL(int)
parse_switches(j_compress_ptr cinfo, int argc, char **argv,
               int last_file_arg_seen, boolean for_real)
/* Parse optional switches.
 * Returns argv[] index of first file-name argument (== argc if none).
 * Any file names with indexes <= last_file_arg_seen are ignored;
 * they have presumably been processed in a previous iteration.
 * (Pass 0 for last_file_arg_seen on the first or only iteration.)
 * for_real is FALSE on the first (dummy) pass; we may skip any expensive
 * processing.
 */
{
  int argn;
  char *arg;
  boolean simple_progressive;
  char *scansarg = NULL;        /* saves -scans parm if any */

  /* Set up default JPEG parameters. */
  simple_progressive = FALSE;
  icc_filename = NULL;
  max_scans = 0;
  outfilename = NULL;
  report = FALSE;
  strict = FALSE;
  copyoption = JCOPYOPT_DEFAULT;
  transformoption.transform = JXFORM_NONE;
  transformoption.perfect = FALSE;
  transformoption.trim = FALSE;
  transformoption.force_grayscale = FALSE;
  transformoption.crop = FALSE;
  transformoption.slow_hflip = FALSE;
  cinfo->err->trace_level = 0;

  /* Scan command line options, adjust parameters */

  for (argn = 1; argn < argc; argn++) {
    arg = argv[argn];
    if (*arg != '-') {
      /* Not a switch, must be a file name argument */
      if (argn <= last_file_arg_seen) {
        outfilename = NULL;     /* -outfile applies to just one input file */
        continue;               /* ignore this name if previously processed */
      }
      break;                    /* else done parsing switches */
    }
    arg++;                      /* advance past switch marker character */

    if (keymatch(arg, "arithmetic", 1)) {
      /* Use arithmetic coding. */
#ifdef C_ARITH_CODING_SUPPORTED
      cinfo->arith_code = TRUE;
#else
      fprintf(stderr, "%s: sorry, arithmetic coding not supported\n",
              progname);
      exit(EXIT_FAILURE);
#endif

    } else if (keymatch(arg, "copy", 2)) {
      /* Select which extra markers to copy. */
      if (++argn >= argc)       /* advance to next argument */
        usage();
      if (keymatch(argv[argn], "none", 1)) {
        copyoption = JCOPYOPT_NONE;
      } else if (keymatch(argv[argn], "comments", 1)) {
        copyoption = JCOPYOPT_COMMENTS;
      } else if (keymatch(argv[argn], "icc", 1)) {
        copyoption = JCOPYOPT_ICC;
      } else if (keymatch(argv[argn], "all", 1)) {
        copyoption = JCOPYOPT_ALL;
      } else
        usage();

    } else if (keymatch(arg, "crop", 2)) {
      /* Perform lossless cropping. */
#if TRANSFORMS_SUPPORTED
      if (++argn >= argc)       /* advance to next argument */
        usage();
      if (transformoption.crop /* reject multiple crop/drop/wipe requests */ ||
          !jtransform_parse_crop_spec(&transformoption, argv[argn])) {
        fprintf(stderr, "%s: bogus -crop argument '%s'\n",
                progname, argv[argn]);
        exit(EXIT_FAILURE);
      }
#else
      select_transform(JXFORM_NONE);    /* force an error */
#endif

    } else if (keymatch(arg, "drop", 2)) {
#if TRANSFORMS_SUPPORTED
      if (++argn >= argc)       /* advance to next argument */
        usage();
      if (transformoption.crop /* reject multiple crop/drop/wipe requests */ ||
          !jtransform_parse_crop_spec(&transformoption, argv[argn]) ||
          transformoption.crop_width_set != JCROP_UNSET ||
          transformoption.crop_height_set != JCROP_UNSET) {
        fprintf(stderr, "%s: bogus -drop argument '%s'\n",
                progname, argv[argn]);
        exit(EXIT_FAILURE);
      }
      if (++argn >= argc)       /* advance to next argument */
        usage();
      dropfilename = argv[argn];
      select_transform(JXFORM_DROP);
#else
      select_transform(JXFORM_NONE);    /* force an error */
#endif

    } else if (keymatch(arg, "debug", 1) || keymatch(arg, "verbose", 1)) {
      /* Enable debug printouts. */
      /* On first -d, print version identification */
      static boolean printed_version = FALSE;

      if (!printed_version) {
        fprintf(stderr, "%s version %s (build %s)\n",
                PACKAGE_NAME, VERSION, BUILD);
        fprintf(stderr, JCOPYRIGHT1);
        fprintf(stderr, JCOPYRIGHT2 "\n");
        fprintf(stderr, "Emulating The Independent JPEG Group's software, version %s\n\n",
                JVERSION);
        printed_version = TRUE;
      }
      cinfo->err->trace_level++;

    } else if (keymatch(arg, "version", 4)) {
      fprintf(stderr, "%s version %s (build %s)\n",
              PACKAGE_NAME, VERSION, BUILD);
      exit(EXIT_SUCCESS);

    } else if (keymatch(arg, "flip", 1)) {
      /* Mirror left-right or top-bottom. */
      if (++argn >= argc)       /* advance to next argument */
        usage();
      if (keymatch(argv[argn], "horizontal", 1))
        select_transform(JXFORM_FLIP_H);
      else if (keymatch(argv[argn], "vertical", 1))
        select_transform(JXFORM_FLIP_V);
      else
        usage();

    } else if (keymatch(arg, "grayscale", 1) ||
               keymatch(arg, "greyscale", 1)) {
      /* Force to grayscale. */
#if TRANSFORMS_SUPPORTED
      transformoption.force_grayscale = TRUE;
#else
      select_transform(JXFORM_NONE);    /* force an error */
#endif

    } else if (keymatch(arg, "icc", 1)) {
      /* Set ICC filename. */
      if (++argn >= argc)       /* advance to next argument */
        usage();
      icc_filename = argv[argn];

    } else if (keymatch(arg, "maxmemory", 3)) {
      /* Maximum memory in Kb (or Mb with 'm'). */
      long lval;
      char ch = 'x';

      if (++argn >= argc)       /* advance to next argument */
        usage();
      if (sscanf(argv[argn], "%ld%c", &lval, &ch) < 1)
        usage();
      if (ch == 'm' || ch == 'M')
        lval *= 1000L;
      cinfo->mem->max_memory_to_use = lval * 1000L;

    } else if (keymatch(arg, "maxscans", 4)) {
      if (++argn >= argc)       /* advance to next argument */
        usage();
      if (sscanf(argv[argn], "%u", &max_scans) != 1)
        usage();

    } else if (keymatch(arg, "optimize", 1) || keymatch(arg, "optimise", 1)) {
      /* Enable entropy parm optimization. */
#ifdef ENTROPY_OPT_SUPPORTED
      cinfo->optimize_coding = TRUE;
#else
      fprintf(stderr, "%s: sorry, entropy optimization was not compiled\n",
              progname);
      exit(EXIT_FAILURE);
#endif

    } else if (keymatch(arg, "outfile", 4)) {
      /* Set output file name. */
      if (++argn >= argc)       /* advance to next argument */
        usage();
      outfilename = argv[argn]; /* save it away for later use */

    } else if (keymatch(arg, "perfect", 2)) {
      /* Fail if there is any partial edge MCUs that the transform can't
       * handle. */
      transformoption.perfect = TRUE;

    } else if (keymatch(arg, "progressive", 2)) {
      /* Select simple progressive mode. */
#ifdef C_PROGRESSIVE_SUPPORTED
      simple_progressive = TRUE;
      /* We must postpone execution until num_components is known. */
#else
      fprintf(stderr, "%s: sorry, progressive output was not compiled\n",
              progname);
      exit(EXIT_FAILURE);
#endif

    } else if (keymatch(arg, "report", 3)) {
      report = TRUE;

    } else if (keymatch(arg, "restart", 1)) {
      /* Restart interval in MCU rows (or in MCUs with 'b'). */
      long lval;
      char ch = 'x';

      if (++argn >= argc)       /* advance to next argument */
        usage();
      if (sscanf(argv[argn], "%ld%c", &lval, &ch) < 1)
        usage();
      if (lval < 0 || lval > 65535L)
        usage();
      if (ch == 'b' || ch == 'B') {
        cinfo->restart_interval = (unsigned int)lval;
        cinfo->restart_in_rows = 0; /* else prior '-restart n' overrides me */
      } else {
        cinfo->restart_in_rows = (int)lval;
        /* restart_interval will be computed during startup */
      }

    } else if (keymatch(arg, "rotate", 2)) {
      /* Rotate 90, 180, or 270 degrees (measured clockwise). */
      if (++argn >= argc)       /* advance to next argument */
        usage();
      if (keymatch(argv[argn], "90", 2))
        select_transform(JXFORM_ROT_90);
      else if (keymatch(argv[argn], "180", 3))
        select_transform(JXFORM_ROT_180);
      else if (keymatch(argv[argn], "270", 3))
        select_transform(JXFORM_ROT_270);
      else
        usage();

    } else if (keymatch(arg, "scans", 1)) {
      /* Set scan script. */
#ifdef C_MULTISCAN_FILES_SUPPORTED
      if (++argn >= argc)       /* advance to next argument */
        usage();
      scansarg = argv[argn];
      /* We must postpone reading the file in case -progressive appears. */
#else
      fprintf(stderr, "%s: sorry, multi-scan output was not compiled\n",
              progname);
      exit(EXIT_FAILURE);
#endif

    } else if (keymatch(arg, "strict", 2)) {
      strict = TRUE;

    } else if (keymatch(arg, "transpose", 1)) {
      /* Transpose (across UL-to-LR axis). */
      select_transform(JXFORM_TRANSPOSE);

    } else if (keymatch(arg, "transverse", 6)) {
      /* Transverse transpose (across UR-to-LL axis). */
      select_transform(JXFORM_TRANSVERSE);

    } else if (keymatch(arg, "trim", 3)) {
      /* Trim off any partial edge MCUs that the transform can't handle. */
      transformoption.trim = TRUE;

    } else if (keymatch(arg, "wipe", 1)) {
#if TRANSFORMS_SUPPORTED
      if (++argn >= argc)       /* advance to next argument */
        usage();
      if (transformoption.crop /* reject multiple crop/drop/wipe requests */ ||
          !jtransform_parse_crop_spec(&transformoption, argv[argn])) {
        fprintf(stderr, "%s: bogus -wipe argument '%s'\n",
                progname, argv[argn]);
        exit(EXIT_FAILURE);
      }
      select_transform(JXFORM_WIPE);
#else
      select_transform(JXFORM_NONE);    /* force an error */
#endif

    } else {
      usage();                  /* bogus switch */
    }
  }

  /* Post-switch-scanning cleanup */

  if (for_real) {

#ifdef C_PROGRESSIVE_SUPPORTED
    if (simple_progressive)     /* process -progressive; -scans can override */
      jpeg_simple_progression(cinfo);
#endif

#ifdef C_MULTISCAN_FILES_SUPPORTED
    if (scansarg != NULL)       /* process -scans if it was present */
      if (!read_scan_script(cinfo, scansarg))
        usage();
#endif
  }

  return argn;                  /* return index of next arg (file name) */
}


METHODDEF(void)
my_emit_message(j_common_ptr cinfo, int msg_level)
{
  if (msg_level < 0) {
    /* Treat warning as fatal */
    cinfo->err->error_exit(cinfo);
  } else {
    if (cinfo->err->trace_level >= msg_level)
      cinfo->err->output_message(cinfo);
  }
}

typedef struct {
    int id;
    int h_samp_fac;
    int w_samp_fac;
    unsigned int block_width;
    unsigned int block_height;
    int dct_scaled_size;
} audrey_component;

typedef struct {
    unsigned int height;
    unsigned int width;
    unsigned int num_components;
    unsigned int colorspace;
    unsigned int bitsofprecision;
    audrey_component components[4];
} audrey_jpeg_header;

__declspec(dllexport)
audrey_jpeg_header
read_jpeg_header(unsigned char* buffer, size_t bufferlen) { // change from int to size_t
    // reformat data into custom struct for py
    audrey_jpeg_header jpgheader;

    // create and intialize struct
    struct jpeg_decompress_struct info;
    struct jpeg_error_mgr jerr;
    info.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&info);

    // process buffer into struct
    jpeg_mem_src(&info, buffer, bufferlen);

    // read in header
    jpeg_read_header(&info, TRUE);

    // copy vals to our struct
    jpgheader.height = info.image_height;
    jpgheader.width = info.image_width;
    jpgheader.num_components = info.num_components;
    jpgheader.colorspace = info.jpeg_color_space;
    jpgheader.bitsofprecision = info.data_precision;

    for (int i = 0; i < info.num_components; i++) {
        audrey_component component;
        component.id = info.cur_comp_info[i]->component_id;
        component.w_samp_fac = info.cur_comp_info[i]->h_samp_factor; // width == horizontal
        component.h_samp_fac = info.cur_comp_info[i]->v_samp_factor; // height == vertical
        component.block_height = info.cur_comp_info[i]->height_in_blocks;
        component.block_width = info.cur_comp_info[i]->width_in_blocks;
        component.dct_scaled_size = info.cur_comp_info[i]->DCT_scaled_size;

        jpgheader.components[i] = component;
    }

    // dealloc struct
    jpeg_destroy_decompress(&info);

    return jpgheader;
}

typedef struct {
    unsigned char* buf;
    unsigned long buflen;
    unsigned long returncode;
} audrey_jpeg;

__declspec(dllexport)
int
free_jpg(unsigned char *jbuf) {
    free(jbuf);
    return 0;
}

__declspec(dllexport)
audrey_jpeg
main(
        int argc, char **argv,
        unsigned char *jsrcbuf, size_t jsrcbuflen,
        unsigned char *jdropbuf, size_t jdropbuflen,
        boolean round_up_partial_mcus)
{
    struct jpeg_decompress_struct srcinfo, dropinfo;
    struct jpeg_compress_struct dstinfo;
    struct jpeg_error_mgr jsrcerr, jdroperr, jdsterr;
    struct cdjpeg_progress_mgr src_progress, dst_progress;
    jvirt_barray_ptr *src_coef_arrays;
    jvirt_barray_ptr *dst_coef_arrays;

    // initialize output struct
    audrey_jpeg dstjpg;
    dstjpg.buf = NULL;
    dstjpg.buflen = 0;

    // REMOVE ME
    /*
    printf("args:\n");
    for (int a = 0;a < argc;a++) {
        printf("%s\n", argv[a]);
    }
    */

    progname = "jpegtran";

    /* Initialize the JPEG decompression object with default error handling. */
    srcinfo.err = jpeg_std_error(&jsrcerr);
    dropinfo.err = jpeg_std_error(&jdroperr);
    jpeg_create_decompress(&srcinfo);

    /* Initialize the JPEG compression object with default error handling. */
    dstinfo.err = jpeg_std_error(&jdsterr);
    jpeg_create_compress(&dstinfo);

    /* Scan command line to find file names.
    * It is convenient to use just one switch-parsing routine, but the switch
    * values read here are mostly ignored; we will rescan the switches after
    * opening the input file.  Also note that most of the switches affect the
    * destination JPEG object, so we parse into that and then copy over what
    * needs to affect the source too.
    */

    parse_switches(&dstinfo, argc, argv, 0, FALSE);
    jsrcerr.trace_level = jdsterr.trace_level;
    srcinfo.mem->max_memory_to_use = dstinfo.mem->max_memory_to_use;

    if (strict)
        jsrcerr.emit_message = my_emit_message;

    if (report) {
        start_progress_monitor((j_common_ptr)&dstinfo, &dst_progress);
        dst_progress.report = report;
    }
    if (report || max_scans != 0) {
        start_progress_monitor((j_common_ptr)&srcinfo, &src_progress);
        src_progress.report = report;
        src_progress.max_scans = max_scans;
    }

    /* Open the drop file. */
    if (jdropbuflen > 0) {
        jpeg_create_decompress(&dropinfo);
        jpeg_mem_src(&dropinfo, jdropbuf, jdropbuflen);
    }

    /* Specify data source for decompression */
    jpeg_mem_src(&srcinfo, jsrcbuf, jsrcbuflen);

    /* Enable saving of extra markers that we want to copy */
    jcopy_markers_setup(&srcinfo, copyoption);

    /* Read file header */
    (void)jpeg_read_header(&srcinfo, TRUE);

    // this rounding equation is shit lol
    if (round_up_partial_mcus) {
        srcinfo.image_height = ((srcinfo.image_height + 15) / 16) * 16;
        srcinfo.image_width = ((srcinfo.image_width + 15) / 16) * 16;
    }

    if (jdropbuflen > 0) {
        (void)jpeg_read_header(&dropinfo, TRUE);

        if (round_up_partial_mcus) {
            dropinfo.image_height = ((dropinfo.image_height + 15) / 16) * 16;
            dropinfo.image_width = ((dropinfo.image_width + 15) / 16) * 16;
        }

        transformoption.crop_width = dropinfo.image_width;
        transformoption.crop_width_set = JCROP_POS;
        transformoption.crop_height = dropinfo.image_height;
        transformoption.crop_height_set = JCROP_POS;
        transformoption.drop_ptr = &dropinfo;
    }

    /* Any space needed by a transform option must be requested before
    * jpeg_read_coefficients so that memory allocation will be done right.
    */
    /* Fail right away if -perfect is given and transformation is not perfect.
    */
    if (!jtransform_request_workspace(&srcinfo, &transformoption)) {
        fprintf(stderr, "%s: transformation is not perfect\n", progname);
        exit(EXIT_FAILURE);
    }

    /* Read source file as DCT coefficients */
    src_coef_arrays = jpeg_read_coefficients(&srcinfo);

    if (jdropbuflen > 0) {
        transformoption.drop_coef_arrays = jpeg_read_coefficients(&dropinfo);
    }

    /* Initialize destination compression parameters from source values */
    jpeg_copy_critical_parameters(&srcinfo, &dstinfo);

    /* Adjust destination parameters if required by transform options;
    * also find out which set of coefficient arrays will hold the output.
    */
    dst_coef_arrays = jtransform_adjust_parameters(&srcinfo, &dstinfo,
                                                    src_coef_arrays,
                                                    &transformoption);

    /* Close input file, if we opened it.
    * Note: we assume that jpeg_read_coefficients consumed all input
    * until JPEG_REACHED_EOI, and that jpeg_finish_decompress will
    * only consume more while (!cinfo->inputctl->eoi_reached).
    * We cannot call jpeg_finish_decompress here since we still need the
    * virtual arrays allocated from the source object for processing.
    */

    /* Adjust default compression parameters by re-parsing the options */
    parse_switches(&dstinfo, argc, argv, 0, TRUE);

    /* Specify data destination for compression */
    jpeg_mem_dest(&dstinfo, &dstjpg.buf, &dstjpg.buflen);
    
    /* Start compressor (note no image data is actually written here) */
    jpeg_write_coefficients(&dstinfo, dst_coef_arrays);

    /* Copy to the output file any extra markers that we want to preserve */
    jcopy_markers_execute(&srcinfo, &dstinfo, copyoption);

    /* Execute image transformation, if any */
    /*
    printf("docrop %i\n", transformoption.crop);
    printf("crop H%i W%i\n", transformoption.output_height, transformoption.output_width);
    */
    jtransform_execute_transformation(&srcinfo, &dstinfo, src_coef_arrays,
                                    &transformoption);

    /* Finish compression and release memory */
    jpeg_finish_compress(&dstinfo);
    jpeg_destroy_compress(&dstinfo);
    (void)jpeg_finish_decompress(&srcinfo);
    jpeg_destroy_decompress(&srcinfo);
    if (jdropbuflen > 0) {
        (void)jpeg_finish_decompress(&dropinfo);
        jpeg_destroy_decompress(&dropinfo);
    }

    if (report)
        end_progress_monitor((j_common_ptr)&dstinfo);
    if (report || max_scans != 0)
        end_progress_monitor((j_common_ptr)&srcinfo);

    /* All done. */
    dstjpg.returncode = (jsrcerr.num_warnings + jdroperr.num_warnings + jdsterr.num_warnings ? EXIT_WARNING : EXIT_SUCCESS);

    return dstjpg;
}

/*
#include <stdio.h>
#include <stdlib.h>
void main() {
    FILE* fp;
    int bufsize;

    // open infile
    fp = fopen(
            "G:\\2023.08.20 realtime sc meta scraper\\scrape_images\\libjpeg-turbo api\\_in.jpg", READ_BINARY);

    // seek to end and get len
    fseek(fp, 0, SEEK_END);
    bufsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    
    unsigned char* bufptr;
    bufptr = malloc(bufsize);
    fread(bufptr, 1, bufsize, fp);

    audrey_jpeg_header balls2 = read_jpeg_header(bufptr, bufsize);

    free(bufptr);
    exit(0);

    unsigned char* balls[] = { "jpegtran", "-perfect", "-copy", "all", "-crop", "104x104", NULL};
    audrey_jpeg penis;

    for (int i = 0; i < 1000; i++) {
        penis = main2(
            6,
            balls,
            bufptr,
            bufsize,
            0,
            0,
            FALSE
        );

        free(penis.buf);

    }

    printf("%i\n", penis.returncode);
    printf("%i\n", penis.buflen);

    free(bufptr);

    exit(0);
}
*/
