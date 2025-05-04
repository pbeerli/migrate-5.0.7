#ifndef ASSIGNMENT_H
#define ASSIGNMENT_H
// assignment of individuals to populations
// works only for a per locus base and then may be summarized
// for all loci, this is similar to allowing for admixture
// each locus is treated indpendently.
// this needs numpop assuming that we know the number of populations
// to simplify the issue with the migration rate matrix
extern void remove_node_assigndb(world_fmt *world,node *p);
extern void fill_world_unassigned(world_fmt *world);
extern void empty_world_unassigned(world_fmt *world);
extern void reset_all_assigned_nodes(world_fmt *world);
extern void reassign_individual(node * node1, long newpop, long otherpop);
extern void report_unassigned(FILE *file, world_fmt *world);
extern void  set_unassigned(node * p, world_fmt * world);
extern void swap_unassigned_nodecollection(world_fmt *world1, world_fmt *world2);
extern void record_assignment(long locus, world_fmt * world);
extern void copy_assignment(world_fmt *world1, world_fmt *world2);
extern void set_unassigned(node * p, world_fmt * world);
extern void update_assignment(world_fmt *world);
extern long chooseUnassigned (proposal_fmt * proposal);
extern void get_assignments (world_fmt * world, option_fmt * options);
extern void swap_unassigned_nodecollection(world_fmt *world1, world_fmt *world2);
extern long find_in_unassignedDB(char *nayme, unassigned_fmt **db);
extern long assign_bypopfreq(world_fmt * world);
extern long assign_bypastfreq(node *thenode, long numpop);
#endif
