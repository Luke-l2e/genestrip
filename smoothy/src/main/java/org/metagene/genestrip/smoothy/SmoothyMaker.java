package org.metagene.genestrip.smoothy;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.smoothy.goals.KMerGenomePathCountsGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.TaxNodeSelection;
import org.metagene.genestrip.tax.TaxTree;

import java.io.File;
import java.util.Map;
import java.util.Set;

public class SmoothyMaker<P extends SmoothyProject> extends GSMaker<P> {
    /**
     * Creates a smoothy maker for the given project.
     *
     * @param project the project to create goals for
     */
    public SmoothyMaker(P project) {
        super(project);
    }

    /**
     * Registers the standard Genestrip goals and additionally the smoothy-specific ones.
     */
    @Override
    protected void createGoals() {
        super.createGoals();

        P project = getProject();


        ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal = (ObjectGoal<Set<RefSeqCategory>, P>) getGoal(GSGoalKey.CATEGORIES);
        ObjectGoal<TaxNodeSelection, P> taxNodesGoal = (ObjectGoal<TaxNodeSelection, P>) getGoal(GSGoalKey.TAXNODES);
        ObjectGoal<TaxTree, P> taxTreeGoal = (ObjectGoal<TaxTree, P>) getGoal(GSGoalKey.TAXTREE);
        RefSeqFnaFilesDownloadGoal fnaFilesGoal = (RefSeqFnaFilesDownloadGoal) getGoal(GSGoalKey.REFSEQFNA);
        ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal = (ObjectGoal<Map<File, TaxTree.TaxIdNode>, P>) getGoal(GSGoalKey.ADD_FASTAS);
        ObjectGoal<AccessionMap, P> accessionMapGoal = (ObjectGoal<AccessionMap, P>) getGoal(GSGoalKey.ACCMAP);
        ObjectGoal<Database, P> dbGoal = (ObjectGoal<Database, P>) getGoal(GSGoalKey.LOAD_DB);

        KMerGenomePathCountsGoal<P> kmerGenomePathCountsGoal = new KMerGenomePathCountsGoal<>(project,
                SmoothyGoalKey.KMER_GENOME_PATH_COUNTS, getExecutionContext(project), categoriesGoal, taxNodesGoal,
                fnaFilesGoal, additionalGoal, accessionMapGoal, dbGoal, taxTreeGoal);
        registerGoal(kmerGenomePathCountsGoal);
    }
}
