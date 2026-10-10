package org.metagene.genestrip.smoothy;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FinerTreeMaker;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.smoothy.goals.KMerGenomePathCountsCSVGoal;
import org.metagene.genestrip.smoothy.goals.KMerGenomePathCountsGoal;
import org.metagene.genestrip.smoothy.goals.KMerGenomePathCountsSERGoal;
import org.metagene.genestrip.smoothy.goals.LoadKMerGenomePathCountsGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.TaxNodeSelection;
import org.metagene.genestrip.tax.TaxTree;

import java.io.File;
import java.util.Map;
import java.util.Set;

public class SmoothyMaker<P extends SmoothyProject> extends FinerTreeMaker<P> {
    /**
     * Creates a smoothy maker for the given project.
     *
     * @param project the project to create goals for
     */
    public SmoothyMaker(P project) {
        super(project);
    }

    /**
     * Registers the standard Genestrip and finer-tree goals and additionally the smoothy-specific ones.
     * <p>
     * The path counts are registered twice, as the finer tree does with its quality goals: once against
     * Genestrip's database and once against the FT database, whose refined nodes give a genome's path
     * more steps. Each variant has goal names of its own, so that both sets of files can sit side by side
     * in the project's {@code db} folder.
     */
    @Override
    protected void createGoals() {
        super.createGoals();

        registerPathCountsGoals((ObjectGoal<Database, P>) getGoal(GSGoalKey.LOAD_DB),
                (FileGoal<P>) getGoal(GSGoalKey.DB), SmoothyGoalKey.KMER_GENOME_PATH_COUNTS,
                SmoothyGoalKey.KMER_GENOME_PATH_COUNTS_CSV, SmoothyGoalKey.KMER_GENOME_PATH_COUNTS_SER,
                SmoothyGoalKey.LOAD_KMER_GENOME_PATH_COUNTS);
        registerPathCountsGoals((ObjectGoal<Database, P>) getGoal(FTGoalKey.LOAD_FTDB),
                (FileGoal<P>) getGoal(FTGoalKey.FTDB), SmoothyGoalKey.FT_KMER_GENOME_PATH_COUNTS,
                SmoothyGoalKey.FT_KMER_GENOME_PATH_COUNTS_CSV, SmoothyGoalKey.FT_KMER_GENOME_PATH_COUNTS_SER,
                SmoothyGoalKey.LOAD_FT_KMER_GENOME_PATH_COUNTS);
    }

    /**
     * Registers the counting, CSV, SER and load goals of the path counts against one database.
     *
     * @param dbGoal      the goal supplying the loaded database
     * @param storeDBGoal the goal storing that database, whose file the loaded counts are checked against
     * @param countsKey   the key of the counting goal
     * @param csvKey      the key of the CSV goal
     * @param serKey      the key of the SER goal
     * @param loadKey     the key of the load goal
     */
    protected void registerPathCountsGoals(ObjectGoal<Database, P> dbGoal, FileGoal<P> storeDBGoal,
                                           SmoothyGoalKey countsKey, SmoothyGoalKey csvKey,
                                           SmoothyGoalKey serKey, SmoothyGoalKey loadKey) {
        P project = getProject();

        ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal = (ObjectGoal<Set<RefSeqCategory>, P>) getGoal(GSGoalKey.CATEGORIES);
        ObjectGoal<TaxNodeSelection, P> taxNodesGoal = (ObjectGoal<TaxNodeSelection, P>) getGoal(GSGoalKey.TAXNODES);
        ObjectGoal<TaxTree, P> taxTreeGoal = (ObjectGoal<TaxTree, P>) getGoal(GSGoalKey.TAXTREE);
        RefSeqFnaFilesDownloadGoal fnaFilesGoal = (RefSeqFnaFilesDownloadGoal) getGoal(GSGoalKey.REFSEQFNA);
        ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal = (ObjectGoal<Map<File, TaxTree.TaxIdNode>, P>) getGoal(GSGoalKey.ADD_FASTAS);
        ObjectGoal<AccessionMap, P> accessionMapGoal = (ObjectGoal<AccessionMap, P>) getGoal(GSGoalKey.ACCMAP);

        KMerGenomePathCountsGoal<P> countsGoal = new KMerGenomePathCountsGoal<>(project, countsKey,
                getExecutionContext(project), categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal,
                accessionMapGoal, dbGoal, taxTreeGoal);
        registerGoal(countsGoal);

        KMerGenomePathCountsCSVGoal<P> csvGoal = new KMerGenomePathCountsCSVGoal<>(project, csvKey, countsGoal,
                dbGoal, getGoal(GSGoalKey.SETUP));
        registerGoal(csvGoal);

        KMerGenomePathCountsSERGoal<P> serGoal = new KMerGenomePathCountsSERGoal<>(project, serKey, countsGoal,
                dbGoal, getGoal(GSGoalKey.SETUP));
        registerGoal(serGoal);

        LoadKMerGenomePathCountsGoal<P> loadGoal = new LoadKMerGenomePathCountsGoal<>(project, loadKey, countsGoal,
                serGoal, storeDBGoal);
        registerGoal(loadGoal);
    }
}
