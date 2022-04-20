/*
 * Copyright© 2022 SemiMod UG. All rights reserved.
 */

#include "ngspice/iferrmsg.h"
#include "ngspice/memory.h"
#include "ngspice/ngspice.h"
#include "ngspice/typedefs.h"

#include "osdidefs.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

/*
 * Handles any errors raised by the setup_instance and setup_model functions
 */
static int handle_init_info(OsdiInitInfo info, const OsdiDescriptor *descr) {
  if (info.flags & (EVAL_RET_FLAG_FATAL | EVAL_RET_FLAG_FINISH)) {
    return (E_PANIC);
  }
  if (info.num_errors == 0) {
    return (OK);
  }

  for (uint32_t i = 0; i < info.num_errors; i++) {
    OsdiInitError *err = &info.errors[i];
    switch (err->code) {
    case INIT_ERR_OUT_OF_BOUNDS: {
      char *param = descr->param_opvar[err->payload.parameter_id].name[0];
      printf("Parameter %s is out of bounds!", param);
      break;
    }
    default:
      printf("Unkown OSDO init error code %d!", err->code);
    }
  }
  free(info.errors);
  errMsg = "Errors occurred during initalization";
  return (E_PRIVATE);
}

/*
 * The OSDI instance data contains the `node_mapping` array.
 * Here an index is stored for each node. This function initalizes this array
 * with its indecies {0, 1, 2, 3, .., n}.
 * The node collapsing information generated by setup_instance is used to
 * replace these initial indecies with those that a node is collapsed into.
 * For example collapsing nodes i and j sets node_mapping[i] = j.
 *
 * Terminals can never be collapsed in ngspice because they are allocated by
 * ngspice instead of OSDI. Therefore any node collapsing that involves nodes
 * `i < connected_terminals` is ignored.
 *
 * @param const OsdiDescriptor *descr The OSDI descriptor
 * @param void *inst The instance data connected_terminals
 * @param uint32_t connected_terminals The number of terminals that are not
 * internal nodes.
 *
 * @returns The number of nodes required after collapsing.
 * */
static uint32_t collapse_nodes(const OsdiDescriptor *descr, void *inst,
                               uint32_t connected_terminals) {
  /* access data inside instance */
  uint32_t *node_mapping =
      (uint32_t *)(((char *)inst) + descr->node_mapping_offset);
  bool *is_collapsible =
      (bool *)(((char *)inst) + descr->is_collapsible_offset);

  /* without collapsing just return the total number of nodes */
  uint32_t num_nodes = descr->num_nodes;

  /*  populate nodes with themselves*/
  for (uint32_t i = 0; i < descr->num_nodes; i++) {
    node_mapping[i] = i;
  }

  for (uint32_t i = 0; i < descr->num_collapsible; i++) {
    /* check if the collapse hint (V(x,y) <+ 0) was executed */
    if (!is_collapsible[i]) {
      continue;
    }

    uint32_t from = descr->collapsible[i].node_1;
    uint32_t to = descr->collapsible[i].node_2;

    /* terminal scan not be collapsed because these are created by the simulator
     */
    if (node_mapping[from] < connected_terminals &&
        (to == descr->num_nodes || node_mapping[to] < connected_terminals ||
         node_mapping[to] == descr->num_nodes)) {
      continue;
    }
    /* ensure that from is always the smaller node */
    if (to != descr->num_nodes && node_mapping[from] < node_mapping[to]) {
      uint32_t temp = from;
      from = to;
      to = temp;
    }

    from = node_mapping[from];
    if (node_mapping[to] != descr->num_nodes) {
      to = node_mapping[to];
    }

    /* replace nodes mapped to from with to and reduce the number of nodes */
    for (uint32_t j = 0; j < descr->num_collapsible; j++) {
      if (node_mapping[j] == from) {
        node_mapping[j] = to;
      } else if (node_mapping[j] > from) {
        node_mapping[j] -= 1;
      }
    }
    num_nodes -= 1;
  }
  return num_nodes;
}

/* replace node mapping local to the current instance (created by
 * collapse_nodes) with global node indicies allocated with CKTmkVolt */
static void write_node_mapping(const OsdiDescriptor *descr, void *inst,
                               uint32_t *nodes) {
  uint32_t *node_mapping =
      (uint32_t *)(((char *)inst) + descr->node_mapping_offset);
  for (uint32_t i = 0; i < descr->num_nodes; i++) {
    if (node_mapping[i] == descr->num_nodes) {
      /* gnd node */
      node_mapping[i] = 0;
    } else {
      node_mapping[i] = nodes[node_mapping[i]];
    }
  }
}

static int init_matrix(SMPmatrix *matrix, const OsdiDescriptor *descr,
                       void *inst) {
  uint32_t *node_mapping =
      (uint32_t *)(((char *)inst) + descr->node_mapping_offset);

  double **jacobian_ptr_resist =
      (double **)(((char *)inst) + descr->jacobian_ptr_resist_offset);

  double **jacobian_ptr_react =
      (double **)(((char *)inst) + descr->jacobian_ptr_react_offset);

  for (uint32_t i = 0; i < descr->num_jacobian_entries; i++) {
    uint32_t equation = descr->jacobian_entries[i].node_1;
    uint32_t unkown = descr->jacobian_entries[i].node_2;
    equation = node_mapping[equation];
    unkown = node_mapping[unkown];
    double *ptr = SMPmakeElt(matrix, (int)equation, (int)unkown);

    if (ptr == NULL) {
      return (E_NOMEM);
    }
    jacobian_ptr_resist[i] = ptr;
    // complex number for ac analysis
    jacobian_ptr_react[i] = ptr + 1;
  }
  return (OK);
}

int OSDIsetup(SMPmatrix *matrix, GENmodel *inModel, CKTcircuit *ckt,
              int *states) {
  OsdiInitInfo init_info;
  OsdiNgspiceHandle handle;
  GENmodel *gen_model;
  int res;
  int error;
  CKTnode *tmp;
  GENinstance *gen_inst;
  int err;

  OsdiRegistryEntry entry = osdi_reg_entry_model(inModel);
  const OsdiDescriptor *descr = entry.descriptor;

  /* setup a temporary buffer */
  uint32_t *node_ids = TMALLOC(uint32_t, descr->num_nodes);

  /* determine the number of states required by each instance */
  int num_states = 0;
  for (uint32_t i = 0; i < descr->num_nodes; i++) {
    if (descr->nodes[i].is_reactive) {
      num_states += 2;
    }
  }

  for (gen_model = inModel; gen_model; gen_model = gen_model->GENnextModel) {
    void *model = osdi_model_data(gen_model);

    /* setup model parameter (setup_model)*/
    handle = (OsdiNgspiceHandle){.kind = 1, .name = gen_model->GENmodName};
    init_info = descr->setup_model((void *)&handle, model);
    res = handle_init_info(init_info, descr);
    if (res) {
      errRtn = "OSDI setup_model";
      continue;
    }

    for (gen_inst = gen_model->GENinstances; gen_inst;
         gen_inst = gen_inst->GENnextInstance) {
      void *inst = osdi_instance_data(entry, gen_inst);

      /* special handling for temperature parameters */
      double temp = ckt->CKTtemp;
      OsdiExtraInstData *extra_inst_data =
          osdi_extra_instance_data(entry, gen_inst);
      if (extra_inst_data->temp_given) {
        temp = extra_inst_data->temp;
      }
      if (extra_inst_data->dt_given) {
        temp += extra_inst_data->dt;
      }

      /* find number of connected ports to allow evaluation of $port_connected
       * and to handle node collapsing correctly later
       * */
      int *terminals = (int *)(gen_inst + 1);
      uint32_t connected_terminals = descr->num_terminals;
      for (uint32_t i = 0; i < descr->num_terminals; i++) {
        if (terminals[i] == -1) {
          connected_terminals = i;
          break;
        }
      }

      /* calculate op independent data, init instance parameters and determine
       which collapsing occurs*/
      handle = (OsdiNgspiceHandle){.kind = 2, .name = gen_inst->GENname};
      init_info = descr->setup_instance((void *)&handle, inst, model, temp,
                                        connected_terminals);
      res = handle_init_info(init_info, descr);
      if (res) {
        errRtn = "OSDI setup_instance";
        continue;
      }

      /* setup the instance nodes */

      uint32_t num_nodes = collapse_nodes(descr, inst, connected_terminals);
      /* copy terminals */
      memcpy(node_ids, gen_inst + 1, sizeof(int) * connected_terminals);
      /* create internal nodes as required */
      for (uint32_t i = connected_terminals; i < num_nodes; i++) {
        // TODO handle currents  correctly
        error = CKTmkVolt(ckt, &tmp, gen_inst->GENname, descr->nodes[i].name);
        if (error)
          return (error);
        node_ids[i] = (uint32_t)tmp->number;
        // TODO nodeset?
      }
      write_node_mapping(descr, inst, node_ids);

      /* now that we have the node mapping we can create the matrix entries */
      err = init_matrix(matrix, descr, inst);
      if (err) {
        return err;
      }

      /* reserve space in the state vector*/
      gen_inst->GENstate = *states;
      *states += num_states;
    }
  }

  free(node_ids);

  return (OK);
}

/* OSDI does not differentiate between setup and temperature update so we just
 * call the setup routines again and assume that  node collapsing (and therefore
 * node mapping) stays the same
 */
extern int OSDItemp(GENmodel *inModel, CKTcircuit *ckt) {
  OsdiInitInfo init_info;
  OsdiNgspiceHandle handle;
  GENmodel *gen_model;
  int res;
  GENinstance *gen_inst;

  OsdiRegistryEntry entry = osdi_reg_entry_model(inModel);
  const OsdiDescriptor *descr = entry.descriptor;

  for (gen_model = inModel; gen_model != NULL;
       gen_model = gen_model->GENnextModel) {
    void *model = osdi_model_data(gen_model);

    handle = (OsdiNgspiceHandle){.kind = 1, .name = gen_model->GENmodName};
    init_info = descr->setup_model((void *)&handle, model);
    res = handle_init_info(init_info, descr);
    if (res) {
      errRtn = "OSDI setup_model (OSDItemp)";
      continue;
    }

    for (gen_inst = gen_model->GENinstances; gen_inst != NULL;
         gen_inst = gen_inst->GENnextInstance) {
      void *inst = osdi_instance_data(entry, gen_inst);

      // special handleing for temperature parameters
      double temp = ckt->CKTtemp;
      OsdiExtraInstData *extra_inst_data =
          osdi_extra_instance_data(entry, gen_inst);
      if (extra_inst_data->temp_given) {
        temp = extra_inst_data->temp;
      }
      if (extra_inst_data->dt_given) {
        temp = extra_inst_data->dt;
      }

      handle = (OsdiNgspiceHandle){.kind = 2, .name = gen_inst->GENname};
      // TODO optional terminals
      init_info = descr->setup_instance((void *)&handle, inst, model, temp,
                                        descr->num_terminals);
      res = handle_init_info(init_info, descr);
      if (res) {
        errRtn = "OSDI setup_instance (OSDItemp)";
        continue;
      }
      // TODO check that there are no changes in node collapse?
    }
  }
  return (OK);
}

/* delete internal nodes
 */
extern int OSDIunsetup(GENmodel *inModel, CKTcircuit *ckt) {
  GENmodel *gen_model;
  GENinstance *gen_inst;
  int num;

  OsdiRegistryEntry entry = osdi_reg_entry_model(inModel);
  const OsdiDescriptor *descr = entry.descriptor;

  for (gen_model = inModel; gen_model != NULL;
       gen_model = gen_model->GENnextModel) {

    for (gen_inst = gen_model->GENinstances; gen_inst != NULL;
         gen_inst = gen_inst->GENnextInstance) {
      void *inst = osdi_instance_data(entry, gen_inst);

      uint32_t *node_mapping =
          (uint32_t *)(((char *)inst) + descr->node_mapping_offset);
      for (uint32_t i = 0; i < descr->num_nodes; i++) {
        num = (int)node_mapping[i];
        // hand coded implementations just know which nodes were collapsed
        // however nodes may be collapsed multiple times so we can't easily use
        // an approach like that instead we delete all nodes
        // Deleting twiche with CKLdltNNum is fine (entry is already removed
        // from the linked list and therefore no action is taken).
        // However CKTdltNNum (rightfully) throws an error when trying to delete
        // an external node. Therefore we need to check for each node that it is
        // an internal node
        if (ckt->prev_CKTlastNode->number &&
            num > ckt->prev_CKTlastNode->number) {
          CKTdltNNum(ckt, num);
        }
      }
    }
  }
  return (OK);
}
