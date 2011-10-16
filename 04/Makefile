# CC=gcc
# CFLAGS=-g -std=c99 -pedantic -Wall

.PHONY: all gv pics grafy kapitoly

all: dh.ps

gv: dh.ps
	(gv $< &)

dh.pdf: dh.tex init_setup.tex fj.sty pictures 
	pdfcslatex $<; bibtex dh; pdfcslatex $<; pdfcslatex $<

dh.ps: dh.dvi
	dvips $<

dh.dvi: dh.tex init_setup.tex fj.sty pictures
	latex -output-format dvi $<; bibtex dh; latex -output-format dvi $<; latex -output-format dvi $<

pictures: pictures/img_5916.jpg pictures/picture-Mitsubishi-Melfa-RV-6S-s_popisky.png pictures/picture-robotic_toolbox_ready_configuration.png
