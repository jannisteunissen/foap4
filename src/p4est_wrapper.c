#include <p4est.h>
#include <p4est_extended.h>
#include <p4est_iterate.h>
#include <p4est_vtk.h>

/* This file provides methods to use p4est from Fortran */

const int FACE_BOUNDARY = -1;
const int FACE_SAME_LEVEL = 0;
const int FACE_COARSE_FINE = 1;

/* Global variables */
p4est_t              *p4est;
p4est_connectivity_t *conn;
sc_MPI_Comm           mpicomm;

int n_faces = 0;

typedef struct bnd_face {
  int face_type;
  int face;                     /* In our geometries, one entry is enough */
  int procs[2];
  int treeid[2];
  int quadid[3];                /* Size 3 for hanging faces */
} bnd_face_t;

bnd_face_t *bnd_face;

void callback_get_faces (p4est_iter_face_info_t * info, void *user_data) {
  p4est_iter_face_side_t *sides = (p4est_iter_face_side_t *) (info->sides.array);

  bnd_face_t *bf = &bnd_face[n_faces];

  bf->treeid[0] = sides[0].treeid;
  bf->quadid[2] = -2;

  if (info->sides.elem_count == 1) {
    bf->face_type = FACE_BOUNDARY;
    bf->face = sides[0].face;

    /* A physical boundary, so there is one full and local face side */
    bf->quadid[0] = sides[0].is.full.quadid;
    bf->procs[0] = sides[0].is.full.is_ghost;
    bf->procs[1] = sides[0].is.full.is_ghost;
    bf->treeid[1] = -2;
    bf->quadid[1] = -2;
  } else if (sides[0].is_hanging | sides[1].is_hanging) {
    /* A coarse-to-fine interface */
    bf->face_type = FACE_COARSE_FINE;

    p4est_iter_face_side_t *face_hanging, *face_full;

    if (sides[0].is_hanging) {
      face_hanging = &sides[0];
      face_full = &sides[1];
    } else {
      face_hanging = &sides[1];
      face_full = &sides[0];
    }

    bf->face = face_full->face;

    /* Store coarse face as first entry */
    bf->quadid[0] = face_full->is.full.quadid;
    bf->procs[0] = face_full->is.full.is_ghost;

    /* Store hanging faces */
    bf->quadid[1] = face_hanging->is.hanging.quadid[0];
    bf->quadid[2] = face_hanging->is.hanging.quadid[1];

    bf->procs[1] = face_hanging->is.hanging.is_ghost[0];

    /* The two hanging faces should be on the same rank */
    P4EST_ASSERT (face_hanging->is.hanging.is_ghost[0] ==
                  face_hanging->is.hanging.is_ghost[1]);
  } else {
    /* Two faces at the same refinement level */
    bf->face_type = FACE_SAME_LEVEL;
    bf->face = sides[0].face;
    bf->treeid[1] = sides[1].treeid;
    bf->quadid[0] = sides[0].is.full.quadid;
    bf->quadid[1] = sides[1].is.full.quadid;

    bf->procs[0] = sides[0].is.full.is_ghost;
    bf->procs[1] = sides[1].is.full.is_ghost;
  }

  n_faces++;
}

/* Initialize MPI and p4est */
void pw_initialize_mpi_and_p4est(MPI_Fint *comm_fortran, int max_blocks) {
  int         argc = 0;
  int         mpiret;

  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Init (&argc, NULL);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  *comm_fortran = MPI_Comm_c2f(mpicomm);

  bnd_face = malloc(max_blocks * 4 * sizeof(bnd_face_t));
}

/* Finalize MPI and p4est */
void pw_finalize_mpi_and_p4est() {
  int mpiret;

  if (p4est != NULL) {
    p4est_destroy (p4est);
  }

  if (conn != NULL) {
    p4est_connectivity_destroy (conn);
  }

  if (bnd_face != NULL) {
    free (bnd_face);
  }

  /* check memory balance and clean up internal registrations */
  sc_finalize ();

  /* release the MPI subsytem */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
}

/* Set brick connectivity */
void pw_set_connectivity_brick(int mi, int ni, int periodic_a, int periodic_b,
                               int min_level, int fill_uniform) {
  conn = p4est_connectivity_new_brick (mi, ni, periodic_a, periodic_b);
  p4est = p4est_new_ext (mpicomm, conn, 0, min_level, fill_uniform,
                         0, NULL, NULL);
}

int pw_get_num_local_quadrants(void) {
  return p4est->local_num_quadrants;
}

/* Get information about local quadrants */
void pw_get_quadrants(int n_quadrants, double *coord, int *level) {
  p4est_topidx_t    tt;
  size_t            zz;
  p4est_tree_t     *tree;
  p4est_quadrant_t *quadrant;
  sc_array_t       *tquadrants;
  int               i_quad;
  double            vxyz[3];

  P4EST_ASSERT (n_quadrants == p4est->local_num_quadrants);
  i_quad = 0;

  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    tquadrants = &tree->quadrants;
    for (zz = 0; zz < tquadrants->elem_count; ++zz) {
      quadrant = p4est_quadrant_array_index (tquadrants, zz);

      p4est_qcoord_to_vertex (conn, tt, quadrant->x, quadrant->y, vxyz);
      coord[2*i_quad] = vxyz[0];
      coord[2*i_quad+1] = vxyz[1];
      level[i_quad] = quadrant->level;
      i_quad++;
    }
  }
}

/* Write mesh to file */
void pw_vtk_write_file(char *fname) {
  p4est_vtk_write_file (p4est, NULL, fname);
}

/* Store information about all faces between quadrants */
void pw_get_all_faces (int *n_faces_arg, bnd_face_t **bnd_face_arg) {
  p4est_ghost_t        *ghost;
  p4est_mesh_t         *mesh;

  const int compute_tree_index = 1;
  const int compute_level_lists = 0;

  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  mesh = p4est_mesh_new_ext(p4est, ghost, compute_tree_index,
                            compute_level_lists, P4EST_CONNECT_FACE);

  /* Determine ghost owners */
  int *ghost_rank = (int *) malloc(sizeof(int) * ghost->ghosts.elem_count);

  for (int i = 0; i < p4est->mpisize; i++) {
    for (int n = ghost->proc_offsets[i]; n < ghost->proc_offsets[i+1]; n++) {
      ghost_rank[n] = i;
    }
  }

  n_faces = 0;

  p4est_iterate (p4est, ghost, NULL, NULL, callback_get_faces, NULL);

  /* Change quadids for ghosts to actual ids */
  for (int n = 0; n < n_faces; n++) {
    for (int i = 0; i < 2; i++) {
      if (bnd_face[n].procs[i]) {
        /* Data from ghost */
        bnd_face[n].procs[i] = ghost_rank[bnd_face[n].quadid[i]];

        int qid = bnd_face[n].quadid[i];
        p4est_quadrant_t *quad = (p4est_quadrant_t *)
          (ghost->ghosts.array + qid * sizeof(p4est_quadrant_t));

        /* This is the quadid cumulative over trees */
        bnd_face[n].quadid[i] = quad->p.piggy3.local_num;
      } else {
        /* Local data */
        bnd_face[n].procs[i] = p4est->mpirank;

        /* Account for tree offset in quadid */
        p4est_tree_t *tree;
        tree = p4est_tree_array_index (p4est->trees, bnd_face[n].treeid[i]);
        bnd_face[n].quadid[i] += tree->quadrants_offset;
      }
    }
  }

  *bnd_face_arg = bnd_face;
  *n_faces_arg = n_faces;

  free (ghost_rank);
  p4est_mesh_destroy (mesh);
  p4est_ghost_destroy (ghost);
}

/* callback to tell p4est which quadrants shall be refined */
static int refine_function (p4est_t *p4est, p4est_topidx_t which_tree,
                          p4est_quadrant_t *quadrant) {
  return (quadrant->p.user_int > 0);
}

/* callback to tell p4est which quadrants shall be coarsened */
static int coarsen_function (p4est_t *p4est, p4est_topidx_t which_tree,
                             p4est_quadrant_t *quadrant[])
{
  int coarsen = 1;

  for (int i = 0; i < P4EST_CHILDREN; ++i) {
    coarsen = (coarsen && quadrant[i]->p.user_int < 0);
  }
  return coarsen;
}

int pw_get_mesh_revision() {
  return p4est->revision;
}

void pw_adjust_refinement(const int n_quadrants, const int *flags,
                          int *has_changed) {
  p4est_topidx_t tt;
  size_t         zz;
  int            i_quad;
  long           old_revision;

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

  /* adapt the new forest non-recursively */
  p4est_refine (p4est, 0, refine_function, NULL);
  p4est_coarsen (p4est, 0, coarsen_function, NULL);

  /* P4EST_CONNECT_FULL is required for filling ghost cells near refinement
     boundaries. */
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  *has_changed = (p4est->revision > old_revision);
}
