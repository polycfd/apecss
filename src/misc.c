// This source file is part of APECSS, an open-source software toolbox
// for the computation of pressure-driven bubble dynamics and acoustic
// emissions in spherical symmetry.
//
// Copyright (C) 2022 The APECSS Developers
//
// The APECSS Developers are listed in the README.md file available in
// the GitHub repository at https://github.com/polycfd/apecss.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "apecss.h"

int apecss_writeonscreen(char* str)
{
  printf("| APECSS | %s\n", str);
  return (0);
}

int apecss_erroronscreen(int num, char* message)
{
  char str[APECSS_STRINGLENGTH_SPRINTF];

  if (num != 0)  // Terminal error. Stop program and exit
  {
    printf("------------------------------------------------------------------------------------ \n");
    sprintf(str, "| APECSS | ERROR: %s", message);
    apecss_writeonscreen(str);
    printf("------------------------------------------------------------------------------------ \n");
    exit(1);
  }
  else
  {
    printf("------------------------------------------------------------------------------------ \n");
    sprintf(str, "| APECSS | WARNING: %s", message);
    apecss_writeonscreen(str);
    printf("------------------------------------------------------------------------------------ \n");
  }

  return (0);
}

int apecss_infoscreen()
{
  printf("| APECSS | Version: v%.1f (%s)\n", APECSS_VERSION_NUM, APECSS_RELEASE_DATE);
#ifdef __OPTIMIZE__
  printf("| APECSS | Compiled in release mode, ");
#else
  printf("| APECSS | Compiled in debug mode, ");
#endif
#if defined(APECSS_PRECISION_LONGDOUBLE)
  printf("with long-double precision \n");
#else
  printf("with double precision \n");
#endif

  return (0);
}

int apecss_readoneoption(FILE* OptionsFile, char* option)
{
  char ch, ch2;
  int line;
  char str[APECSS_STRINGLENGTH_SPRINTF];

  line = 1;
  while ((ch = getc(OptionsFile)) != EOF)
  {
    if (ch != '\n' && ch != '#')
    {
      while (ch == ' ' || ch == '\t' || ch == '\n')
      {
        ch = getc(OptionsFile);
      }

      memcpy(option, &ch, sizeof(char));
      ch2 = getc(OptionsFile);
      ungetc(ch2, OptionsFile);

      if (ch2 != ' ' && ch2 != '\n')
      {
        apecss_linegetskip(str, OptionsFile);
        memcpy(option + 1, str, sizeof(str));
      }
      else
      {
        str[0] = '\0';
        option[1] = '\0';
      }

      return (line);
    }

    if (ch != '\n')
    {
      apecss_lineget(str, OptionsFile);
    }

    line++;
  }

  return (EOF);
}

int apecss_lineget(char* str, FILE* OptionsFile)
{
  char ch;
  int index = 0;

  while ((ch = getc(OptionsFile)) != EOF)
  {
    str[index++] = ch;
    if (index == APECSS_STRINGLENGTH_SPRINTF) apecss_erroronscreen(1, "String in file too long!");

    if (ch == '\n')
    {
      str[index] = '\0';
      return (index - 1);
    }
  }

  return (0);
}

int apecss_linegetskip(char* str, FILE* OptionsFile)
{
  char ch;
  int index = 0;

  while ((ch = getc(OptionsFile)) != EOF)
  {
    if (ch != ' ' && ch != '\n')
    {
      str[index++] = ch;
      if (index == APECSS_STRINGLENGTH_SPRINTF) apecss_erroronscreen(1, "String in file too long!");
    }
    if ((ch == '\n' || ch == ' ') && index > 0)
    {
      str[index] = '\0';
      return (index);
    }
  }

  return (0);
}