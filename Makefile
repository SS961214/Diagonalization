FRC = force_rebuild # If you set FRC=force_rebuild, dependencies lists will be newly generated.
PREFIX ?= $(HOME)
TARGETROOT = $(PREFIX:%/=%)/bin/Diagonalization

SUBDIRS = ./Headers $(shell find ./*/* -name Makefile -or -name basic.mk | grep -v -e "Headers" | sed -e 's/\/Makefile//' -e 's/\/basic.mk//' | sort | uniq)
EXECS   = $(shell find ./*/* -name "*.out")
OBJS    = $(shell find ./*/* -name "*.o")
TOCOPY  = 
SHELLSCRIPTS = $(shell find ./* -name "*.sh")
SHELLSCRIPTS += $(shell find ./* -name "CreateScript*.c")
DIRS = $(dir $(EXECS)) $(dir $(SHELLSCRIPTS))

.PHONY: test all subdirs $(SUBDIRS)

all : subdirs $(SHELLSCRIPTS)
	@ $(RM) -r $(TARGETROOT)
	@for out in $(DIRS:./%=%) ;do mkdir -p $(TARGETROOT)/$$out;done
	@for out in $(EXECS:./%=%);do mv $$out $(TARGETROOT)/$$out;done
	@for out in $(SHELLSCRIPTS:./%=%);do cp $$out $(TARGETROOT)/$$out;done
	@for out in $(TOCOPY)     ;do cp -rf $$out $(TARGETROOT)/$$out;done
	@for objs in $(OBJS);do $(RM) $$objs;done
	@echo
	@COLOR_f=$$(printf '\033[36m'); COLOR_b=$$(printf '\033[m'); echo "$${COLOR_f}(Make completed) Executable files are placed in \"$(TARGETROOT:%/=%)\".$${COLOR_b}"
	@COLOR_f=$$(printf '\033[36m'); COLOR_b=$$(printf '\033[m'); echo "$${COLOR_f}If you want to specify other places, use 'make PREFIX=(target dir).'$${COLOR_b}"
	@echo

subdirs : $(SUBDIRS)

$(SUBDIRS) :
	@echo
	-$(MAKE) -C $@ -f basic.mk Makefile "FRC=${FRC}"
	 $(MAKE) -C $@ -f Makefile

clean :
	@for execs in $(EXECS)  ;do $(RM) $$execs         ;echo "rm -f $$execs"         ;done
	@for objs  in $(OBJS)   ;do $(RM) $$objs          ;echo "rm -f $$objs"          ;done
	@for makes in $(SUBDIRS);do $(RM) $$makes/Makefile;echo "rm -f $$makes/Makefile";done

test :
	for dir in $(SUBDIRS);do echo $$dir;done
	@for out in $(DIRS:./%=%)        ;do mkdir -p $(TARGETROOT)/$$out;done
	@for out in $(EXECS:./%=%)       ;do mv $$out $(TARGETROOT)/$$out;done
	@for out in $(SHELLSCRIPTS:./%=%);do cp $$out $(TARGETROOT)/$$out;done
