LHS_SRC = Complex.lhs StepGen.lhs YeeLattice.lhs

HS_SRC = step_d_gen.hs step_h_gen.hs update_e_from_d_gen.hs update_from_e_gen.hs

SRC = $(LHS_SRC) $(HS_SRC)

EXTRA_DIST = .depend $(SRC)

HC_FLAGS = -fglasgow-exts -package data

noinst_SCRIPTS = step_d_gen step_h_gen update_e_from_d_gen update_from_e_gen

SUFFIXES = .o .hs .lhs .hi

.o.hi:
	@:

.hs.o:
	$(HC) -c -o $@ $(HC_FLAGS) $<

.lhs.o:
	$(HC) -c -o $@ $(HC_FLAGS) $<

STEP_H_GEN_OBJ = step_h_gen.o StepGen.o YeeLattice.o
step_h_gen: $(STEP_H_GEN_OBJ)
	$(HC) -o $@ $(HC_FLAGS) $(STEP_H_GEN_OBJ)

STEP_D_GEN_OBJ = step_d_gen.o StepGen.o YeeLattice.o
step_d_gen: $(STEP_D_GEN_OBJ)
	$(HC) -o $@ $(HC_FLAGS) $(STEP_D_GEN_OBJ)

UPDATE_E_FROM_D_GEN_OBJ = update_e_from_d_gen.o Complex.o StepGen.o YeeLattice.o
update_e_from_d_gen: $(UPDATE_E_FROM_D_GEN_OBJ)
	$(HC) -o $@ $(HC_FLAGS) $(UPDATE_E_FROM_D_GEN_OBJ)

UPDATE_FROM_E_GEN_OBJ = update_from_e_gen.o Complex.o StepGen.o YeeLattice.o
update_from_e_gen: $(UPDATE_FROM_E_GEN_OBJ)
	$(HC) -o $@ $(HC_FLAGS) $(UPDATE_FROM_E_GEN_OBJ)

clean-local:
	rm -f *.o *.hi $(noinst_SCRIPTS)

depend:
	(cd ${srcdir}; $(HC) -M -optdep-f -optdep.depend $(HC_FLAGS) $(SRC))

include .depend
