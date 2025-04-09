#include <p4est.h>
#include <p4est_extended.h>
#include <p4est_iterate.h>
#include <p4est_vtk.h>
#include <p4est_communication.h>

/* This file provides methods to use p4est from Fortran */

const int FACE_BOUNDARY = 0;
const int FACE_SAME_LEVEL = 1;
const int FACE_COARSE_TO_FINE = 2;
const int FACE_FINE_TO_COARSE = 3;

/* Struct to store information about a quadrant face */
typedef struct bnd_face {
  int face_type;                /* What kind of face (same level, ...) */
  int face;                     /* Direction of the face */
  int other_proc;               /* MPI rank that owns quadid[1] */
  int quadid[2];                /* quadid[0] is always local, [1] can be non-local */
  int offset;                   /* Offset for a hanging face */
  int ibuf_recv;                /* Index in receive buffer (not filled here) */
  int ibuf_send;                /* Index in send buffer (not filled here) */
} bnd_face_t;

/* Stores p4est and related data required for this wrapper */
typedef struct pw_state {
  p4est_t              *p4est;
  p4est_connectivity_t *conn;
  p4est_ghost_t        *ghost;
  int                  *ghost_rank;
  sc_MPI_Comm           mpicomm;
  int                   n_faces;
  int                   max_n_faces;
  bnd_face_t           *bnd_face;
} pw_state_t;

/* Initialize MPI and p4est */
void pw_initialize(pw_state_t **pw_ptr, MPI_Fint *comm_fortran,
                   int log_level) {

  int         argc = 0;
  int         mpiret;
  pw_state_t *pw;

  pw = (pw_state_t *) malloc(sizeof(pw_state_t));
  pw->mpicomm = sc_MPI_COMM_WORLD;

  mpiret = sc_MPI_Init (&argc, NULL);
  SC_CHECK_MPI (mpiret);

  sc_init (pw->mpicomm, 1, 1, NULL, log_level);
  p4est_init (NULL, log_level);

  *comm_fortran = MPI_Comm_c2f(pw->mpicomm);

  *pw_ptr = pw;
}

/* Destroy p4est data */
void pw_destroy(pw_state_t *pw) {
  if (pw->p4est != NULL) {
    p4est_destroy (pw->p4est);
  }

  if (pw->conn != NULL) {
    p4est_connectivity_destroy (pw->conn);
  }

  if (pw->bnd_face != NULL) {
    free (pw->bnd_face);
  }
}

/* Finalize MPI and p4est */
void pw_finalize(pw_state_t *pw) {
  int mpiret;

  /* check memory balance and clean up internal registrations */
  sc_finalize ();

  /* release the MPI subsytem */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
}

/* Set brick connectivity */
void pw_set_connectivity_brick(pw_state_t *pw, const int mi, const int ni,
                               const int periodic_a, const int periodic_b,
                               const int min_level, const int fill_uniform,
                               const int max_blocks) {
  pw->conn = p4est_connectivity_new_brick (mi, ni, periodic_a, periodic_b);
  pw->p4est = p4est_new_ext (pw->mpicomm, pw->conn, 0, min_level, fill_uniform,
                         0, NULL, NULL);

  /* Allocate storage for face boundaries */
  pw->max_n_faces = max_blocks * 4 + 10;
  pw->bnd_face = malloc(pw->max_n_faces * sizeof(bnd_face_t));
}

int pw_get_num_local_quadrants(pw_state_t *pw) {
  return pw->p4est->local_num_quadrants;
}

int pw_get_num_global_quadrants(pw_state_t *pw) {
  return pw->p4est->global_num_quadrants;
}

/* Get information about local quadrants */
void pw_get_quadrants(pw_state_t *pw, int n_quadrants,
                      double *coord, int *level) {
  p4est_topidx_t    tt;
  size_t            zz;
  p4est_tree_t     *tree;
  p4est_quadrant_t *quadrant;
  sc_array_t       *tquadrants;
  int               i_quad;
  double            vxyz[3];
  p4est_t           *p4est;

  p4est = pw->p4est;
  P4EST_ASSERT (n_quadrants == p4est->local_num_quadrants);
  i_quad = 0;

  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    tquadrants = &tree->quadrants;
    for (zz = 0; zz < tquadrants->elem_count; ++zz) {
      quadrant = p4est_quadrant_array_index (tquadrants, zz);

      p4est_qcoord_to_vertex (pw->conn, tt, quadrant->x, quadrant->y, vxyz);
      coord[2*i_quad] = vxyz[0];
      coord[2*i_quad+1] = vxyz[1];
      level[i_quad] = quadrant->level;
      i_quad++;
    }
  }
}

/* Write mesh to file */
void pw_vtk_write_file(pw_state_t *pw, char *fname) {
  p4est_vtk_write_file (pw->p4est, NULL, fname);
}

/* Callback function to store a list with all the faces between quadrants */
void callback_get_faces (p4est_iter_face_info_t * info, void *user_data) {
  int                     i, j;
  sc_array_t             *trees;
  int                    *ghost_rank;
  p4est_tree_t           *tree;
  p4est_iter_face_side_t *sides;
  bnd_face_t             *bf;
  int                     ghost_ix;
  p4est_quadrant_t       *quad;
  pw_state_t             *pw;

  pw = (pw_state_t *) user_data;
  trees = (sc_array_t *) pw->p4est->trees;
  ghost_rank = (int *) pw->ghost_rank;
  sides = (p4est_iter_face_side_t *) (info->sides.array);

  if (pw->n_faces > pw->max_n_faces - 10) {
    SC_ABORT("Too many faces in callback_get_faces()");
  }

  if (info->sides.elem_count == 1) {
    /* A physical boundary, so there is just one full (and local) face side */
    bf = &(pw->bnd_face)[pw->n_faces++];
    bf->face_type = FACE_BOUNDARY;
    bf->face = sides[0].face;

    tree = p4est_tree_array_index (trees, sides[0].treeid);
    bf->quadid[0] = sides[0].is.full.quadid + tree->quadrants_offset;
    bf->other_proc = pw->p4est->mpirank;
    bf->quadid[1] = -2;
    bf->offset = -2;
  } else if (sides[0].is_hanging | sides[1].is_hanging) {
    /* A coarse-to-fine interface */
    p4est_iter_face_side_t *face_hanging, *face_full;

    if (sides[0].is_hanging) {
      face_hanging = &sides[0];
      face_full = &sides[1];
    } else {
      face_hanging = &sides[1];
      face_full = &sides[0];
    }

    /* Handle non-ghost hanging faces */
    for (int i = 0; i<2; i++) {
      if (!face_hanging->is.hanging.is_ghost[i]) {
        /* Hanging face is not a ghost, coarse side can be a ghost */
        bf = &(pw->bnd_face)[pw->n_faces++];
        bf->face_type = FACE_FINE_TO_COARSE;
        bf->face = face_hanging->face;

        tree = p4est_tree_array_index (trees, face_hanging->treeid);
        bf->quadid[0] = face_hanging->is.hanging.quadid[i] +
          tree->quadrants_offset;
        bf->offset = i;

        if (face_full->is.full.is_ghost) {
          ghost_ix = face_full->is.full.quadid;
          quad = (p4est_quadrant_t *) (pw->ghost->ghosts.array +
                                       ghost_ix * sizeof(p4est_quadrant_t));
          bf->other_proc = ghost_rank[ghost_ix];
          bf->quadid[1] = quad->p.piggy3.local_num;
        } else {
          bf->other_proc = pw->p4est->mpirank;
          tree = p4est_tree_array_index (trees, face_full->treeid);
          bf->quadid[1] = face_full->is.full.quadid + tree->quadrants_offset;
        }
      }
    }

    /* Handle 'ghost' hanging faces */
    if (!face_full->is.full.is_ghost) {
      for (int i = 0; i<2; i++) {
        if (face_hanging->is.hanging.is_ghost[i]) {
          bf = &(pw->bnd_face)[pw->n_faces++];
          bf->face_type = FACE_COARSE_TO_FINE;
          bf->face = face_full->face;

          tree = p4est_tree_array_index (trees, face_full->treeid);
          bf->quadid[0] = face_full->is.full.quadid + tree->quadrants_offset;
          bf->offset = i;

          ghost_ix = face_hanging->is.hanging.quadid[i];
          quad = (p4est_quadrant_t *) (pw->ghost->ghosts.array +
                                       ghost_ix * sizeof(p4est_quadrant_t));
          bf->other_proc = ghost_rank[ghost_ix];
          bf->quadid[1] = quad->p.piggy3.local_num;
        }
      }
    }
  } else {
    /* Boundary at the same refinement level */
    bf = &(pw->bnd_face)[pw->n_faces++];

    /* Store local side first */
    i = sides[0].is.full.is_ghost; /* 1 or 0 */
    j = !i;

    /* Two faces at the same refinement level */
    bf->face_type = FACE_SAME_LEVEL;
    bf->face      = sides[i].face;

    tree = p4est_tree_array_index (trees, sides[i].treeid);
    bf->quadid[0] = sides[i].is.full.quadid + tree->quadrants_offset;

    if (sides[j].is.full.is_ghost) {
      ghost_ix = sides[j].is.full.quadid;
      quad = (p4est_quadrant_t *) (pw->ghost->ghosts.array +
                                   ghost_ix * sizeof(p4est_quadrant_t));

      bf->other_proc = ghost_rank[ghost_ix];
      bf->quadid[1] = quad->p.piggy3.local_num;

    } else {
      bf->other_proc = pw->p4est->mpirank;
      tree = p4est_tree_array_index (trees, sides[j].treeid);
      bf->quadid[1] = sides[j].is.full.quadid + tree->quadrants_offset;
    }
    bf->offset = -2;
  }
}

/* Store information about all faces between quadrants */
void pw_get_all_faces (pw_state_t *pw, int *n_faces_arg,
                       bnd_face_t **bnd_face_arg) {
  pw->ghost = p4est_ghost_new (pw->p4est, P4EST_CONNECT_FACE);

  /* Determine ghost owners */
  pw->ghost_rank = (int *) malloc(sizeof(int) * pw->ghost->ghosts.elem_count);

  for (int i = 0; i < pw->p4est->mpisize; i++) {
    for (int n = pw->ghost->proc_offsets[i];
         n < pw->ghost->proc_offsets[i+1]; n++) {
      pw->ghost_rank[n] = i;
    }
  }

  pw->n_faces = 0;
  p4est_iterate (pw->p4est, pw->ghost, pw, NULL, callback_get_faces, NULL);

  *bnd_face_arg = pw->bnd_face;
  *n_faces_arg = pw->n_faces;

  free (pw->ghost_rank);
  p4est_ghost_destroy (pw->ghost);
}

/* Callback to tell p4est which quadrants shall be refined */
static int refine_function (p4est_t *p4est, p4est_topidx_t which_tree,
                          p4est_quadrant_t *quadrant) {
  return (quadrant->p.user_int > 0);
}

/* Callback to tell p4est which quadrants shall be coarsened */
static int coarsen_function (p4est_t *p4est, p4est_topidx_t which_tree,
                             p4est_quadrant_t *quadrant[])
{
  int coarsen = 1;

  for (int i = 0; i < P4EST_CHILDREN; ++i) {
    coarsen = (coarsen && quadrant[i]->p.user_int < 0);
  }
  return coarsen;
}

/* Return the mesh revision number, which is incremented when the
   mesh or partition changes */
int pw_get_mesh_revision(pw_state_t *pw) {
  return pw->p4est->revision;
}

/* Get highest local level */
int pw_get_highest_local_level(pw_state_t *pw) {
  p4est_topidx_t  tt;
  p4est_t        *p4est    = pw->p4est;
  int             maxlevel = 0;

  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    p4est_tree_t *tree = p4est_tree_array_index (p4est->trees, tt);
    if (tree->maxlevel > maxlevel) {
      maxlevel = tree->maxlevel;
    }
  }
  return maxlevel;
}

/* Update the refinement non-recursively based on refinement flags */
void pw_adjust_refinement(pw_state_t *pw, const int n_quadrants,
                          const int  *flags, int *has_changed) {
  p4est_topidx_t  tt;
  size_t          zz;
  int             i_quad;
  long            old_revision;
  p4est_t        *p4est = pw->p4est;

  P4EST_ASSERT (n_quadrants == p4est->local_num_quadrants);
  old_revision = p4est->revision;
  i_quad = 0;

  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    p4est_tree_t *tree = p4est_tree_array_index (p4est->trees, tt);
    sc_array_t *tquadrants = &tree->quadrants;
    for (zz = 0; zz < tquadrants->elem_count; ++zz) {
      p4est_quadrant_t *quadrant = p4est_quadrant_array_index (tquadrants, zz);
      quadrant->p.user_int = flags[i_quad];
      i_quad++;
    }
  }

  /* Adapt the new forest non-recursively */
  p4est_refine (p4est, 0, refine_function, NULL);
  p4est_coarsen (p4est, 0, coarsen_function, NULL);

  /* P4EST_CONNECT_FULL is required to ensure ghost cells near refinement
     boundaries are always filled correctly */
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  *has_changed = (p4est->revision > old_revision);
}

/* Re-partition the quadrants over the MPI ranks */
void pw_partition(pw_state_t *pw, int *n_changed, int64_t *gfq_old,
                  int64_t *gfq_new) {
  /* Ensures siblings are at the same MPI rank */
  const int allow_for_coarsening = 1;

  for (int rank = 0; rank < pw->p4est->mpisize+1; rank++) {
    gfq_old[rank] = pw->p4est->global_first_quadrant[rank];
  }

  *n_changed = p4est_partition_ext (pw->p4est, allow_for_coarsening, NULL);

  for (int rank = 0; rank < pw->p4est->mpisize+1; rank++) {
    gfq_new[rank] = pw->p4est->global_first_quadrant[rank];
  }
}

/* Transfer user data after partitioning */
void pw_partition_transfer(pw_state_t *pw, const int64_t *gfq_old,
                           const void *src_data, void *dest_data,
                           const int data_size) {
  const int tag = 0;

  p4est_transfer_fixed (pw->p4est->global_first_quadrant, gfq_old,
                        pw->p4est->mpicomm, tag,
                        dest_data, src_data, data_size);
}
