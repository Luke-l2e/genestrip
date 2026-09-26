/*
 *
 * “Commons Clause” License Condition v1.0
 *
 * The Software is provided to you by the Licensor under the License,
 * as defined below, subject to the following condition.
 *
 * Without limiting other conditions in the License, the grant of rights under the License
 * will not include, and the License does not grant to you, the right to Sell the Software.
 *
 * For purposes of the foregoing, “Sell” means practicing any or all of the rights granted
 * to you under the License to provide to third parties, for a fee or other consideration
 * (including without limitation fees for hosting or consulting/ support services related to
 * the Software), a product or service whose value derives, entirely or substantially, from the
 * functionality of the Software. Any license notice or attribution required by the License
 * must also include this Commons Clause License Condition notice.
 *
 * Software: genestrip
 *
 * License: Apache 2.0
 *
 * Licensor: Daniel Pfeifer (daniel.pfeifer@progotec.de)
 *
 */
package org.metagene.genestrip.smoothy.goals;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.goals.refseq.FastaReaderGoal;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.*;
import org.metagene.genestrip.smoothy.SmoothyProject;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.tax.TaxNodeSelection;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Goal computing, for every genome of the database, how many of its k-mer occurrences the database
 * holds at each node from the genome's own node up to the root.
 * <p>
 * It re-reads every reference genome and, for every k-mer position (subject only to DUST filtering and
 * {@code kMerSampling}, exactly as {@link org.metagene.genestrip.goals.refseq.KMerRankStatsGoal} does),
 * looks the k-mer up directly in the already-built {@link Database}'s {@link KMerStore}. No Bloom
 * filter is involved and nothing is deduplicated: a k-mer occurring twice in a genome is counted twice,
 * which is what distinguishes this walk from the one that fills the database.
 * <p>
 * The result is a {@code Map<String, long[]>} keyed by a genome's tax id (a {@code GENOME}-rank node,
 * requiring the database to have been built with {@code genomeNodes=true}). Each value is an array
 * indexed by position along the genome's own path to the root: index {@code 0} is the genome's own
 * node, and each following index is one ancestor further up. A k-mer's occurrence is tallied at the
 * index of the node the database stores it at, which - by construction of the lowest-common-ancestor
 * update the database was built with - always lies somewhere on that path.
 *
 * @param <P> the project type
 */
public class KMerGenomePathCountsGoal<P extends SmoothyProject> extends FastaReaderGoal<Map<String, long[]>, P> {
    private final ObjectGoal<AccessionMap, P> accessionMapGoal;
    private final ObjectGoal<Database, P> dbGoal;

    // Populated for the duration of doMakeThis() and read by the (multi-threaded) fasta readers.
    private KMerStore<String> store;
    private SmallTaxTree taxTree;
    private Map<String, long[]> map;

    /**
     * Creates the goal, wiring the accession map and the loaded database it looks k-mers up in, in
     * addition to the standard FASTA-reader dependencies.
     *
     * @param project          the project this goal belongs to
     * @param key              the key identifying this goal
     * @param bundle           the execution context supplying the worker threads
     * @param categoriesGoal   the goal supplying the RefSeq categories to read
     * @param taxNodesGoal     the goal supplying the required taxonomy nodes
     * @param fnaFilesGoal     the goal supplying the downloaded RefSeq {@code .fna} files
     * @param additionalGoal   the goal supplying additional FASTA files mapped to their tax node
     * @param accessionMapGoal the goal supplying the accession-to-taxid map
     * @param dbGoal           the goal supplying the loaded database to look k-mers up in
     * @param taxTreeGoal      the goal supplying the taxonomy tree; kept as a direct dependency so that
     *                         the (aggressively cleaned) tree is not freed and rebuilt as a second
     *                         instance while this goal runs, same reason as in {@code KMerRankStatsGoal}
     * @param deps             any further goals this goal depends on
     */
    @SafeVarargs
    public KMerGenomePathCountsGoal(P project, GoalKey key, ExecutionContext bundle,
                                    ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                                    ObjectGoal<TaxNodeSelection, P> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                                    ObjectGoal<Map<File, TaxIdNode>, P> additionalGoal,
                                    ObjectGoal<AccessionMap, P> accessionMapGoal, ObjectGoal<Database, P> dbGoal,
                                    ObjectGoal<TaxTree, P> taxTreeGoal, Goal<P>... deps) {
        super(project, key, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal,
                Goal.append(deps, accessionMapGoal, dbGoal, taxTreeGoal));
        this.accessionMapGoal = accessionMapGoal;
        this.dbGoal = dbGoal;
    }

    @Override
    protected void doMakeThis() {
        try {
            map = new HashMap<>();
            Database database = dbGoal.get();
            store = database.getKmerStore();
            taxTree = database.getTaxTree();
            // One entry per genome, pre-created here (before the multi-threaded readFastas() run) so
            // that the map is never structurally modified while the reader threads work - they only
            // mutate an already-present long[] value, each guarded by its genome node's own monitor.
            for (SmallTaxIdNode node : taxTree) {
                if (Rank.GENOME.equals(node.getRank())) {
                    int pathLength = 0;
                    for (SmallTaxIdNode a = node; a != null; a = a.getParent()) {
                        pathLength++;
                    }
                    map.put(node.getTaxId(), new long[pathLength]);
                }
            }
            if (map.isEmpty() && getLogger().isWarnEnabled()) {
                // Otherwise the goal reads every genome only to hand back an empty map without a word.
                getLogger().warn("The database holds no GENOME-rank nodes - was it built with genomeNodes=true?");
            }
            readFastas();
            set(map);
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            map = null;
            store = null;
            taxTree = null;
        }
    }

    @Override
    protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie contigsPerTaxid) {
        return new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES), taxNodesGoal.get().getSelected(),
                isIncludeRefSeqFna() ? accessionMapGoal.get() : null, intConfigValue(GSConfigKey.KMER_SIZE), intConfigValue(GSConfigKey.MAX_DUST),
                intConfigValue(GSConfigKey.KMER_SAMPLING), booleanConfigValue(GSConfigKey.ASSEMBLY_ACCESSIONS_ONLY),
                contigsPerTaxid, booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES),
                (Rank) configValue(GSConfigKey.FOLD_TAXA_BELOW));
    }

    /**
     * Returns the node the database build filed a genome at, given the node its accession names. A copy
     * of {@code ReworkingStoreFastaReader.foldUp()}, which core keeps package-private: it walks up to
     * the first ancestor of rank {@code below}, never past it, and only if that ancestor is requested.
     * The two must stay in step, or the goal looks for genomes where the build did not put them.
     *
     * @param node     the node the accession resolved to
     * @param below    the rank the build folded to, or {@code null} for no folding
     * @param taxNodes the requested tax nodes, empty for no restriction
     * @return the node the genome was filed at
     */
    static TaxIdNode foldUp(TaxIdNode node, Rank below, Set<TaxIdNode> taxNodes) {
        if (below == null || node == null) {
            return node;
        }
        for (TaxIdNode n = node; n != null; n = n.getParent()) {
            Rank rank = n.getRank();
            if (rank == null || !rank.isComparableTo(below)) {
                continue;
            }
            if (rank == below) {
                return taxNodes == null || taxNodes.isEmpty() || taxNodes.contains(n) ? n : node;
            }
            if (rank.isAbove(below)) {
                // Past the rank asked for without having met it.
                return node;
            }
        }
        return node;
    }

    /**
     * FASTA reader that, for each k-mer of a contig, looks the k-mer up in the database and tallies it
     * against the current contig's genome, at the position along the genome's own path the stored node
     * sits at.
     */
    protected class MyFastaReader extends AbstractStoreFastaReader {
        // Resolved once per contig by infoLine(), not once per k-mer - resolving it is a short walk of
        // the database's compact tree, and a contig's genome does not change while it is being read.
        private SmallTaxIdNode genomeNode;
        private long[] genomeCounts;
        // The rank the database build folded genomes up to, or null - see foldUp().
        // Needed to find the node the build filed a genome under.
        private final Rank foldTaxaBelow;

        /**
         * Creates the reader.
         *
         * @param bufferSize             the read buffer size
         * @param taxNodes               the requested tax nodes
         * @param accessionMap           the accession-to-taxid map
         * @param k                      the k-mer length
         * @param maxDust                the maximum allowed low-complexity (dust) run length
         * @param kMerSampling           the k-mer sampling step size
         * @param assemblyAccessionsOnly whether only genomic accessions are considered, dropping `NG_`, `NT_` and `NW_`
         * @param contigsPerTaxid        the per-taxid contig trie
         * @param enableLowerCaseBases   whether lower-case bases are included
         * @param foldTaxaBelow          the rank the database build folded genomes up to, or {@code null}
         */
        public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k,
                             int maxDust, int kMerSampling, boolean assemblyAccessionsOnly,
                             StringLong2DigitTrie contigsPerTaxid, boolean enableLowerCaseBases,
                             Rank foldTaxaBelow) {
            super(bufferSize, taxNodes, accessionMap, k, maxDust, kMerSampling, assemblyAccessionsOnly,
                    contigsPerTaxid, enableLowerCaseBases);
            this.foldTaxaBelow = foldTaxaBelow;
        }

        /**
         * Resolves and caches this contig's genome node and its counter array, in addition to the
         * standard per-contig resolution the superclass already does.
         */
        @Override
        protected void infoLine() {
            super.infoLine();
            updateGenomeNode();
        }

        /**
         * Resolves {@link #genomeNode} for the current contig: starting from the tax node
         * {@link #infoLine()} already resolved, descends into the database's compact tree to the
         * {@code GENOME}-rank descendant the contig's accession names, mirroring the same genome-key
         * cut {@code ReworkingStoreFastaReader} used when the database was built. Leaves both fields
         * {@code null} when the contig is not included, or the database holds no genome node for it -
         * e.g. because it was not built with {@code genomeNodes=true}.
         */
        private void updateGenomeNode() {
            genomeNode = null;
            genomeCounts = null;
            if (!includeContig || node == null || node.getTaxId() == null) {
                return;
            }
            // The build filed the genome under the folded node, not the one the accession names: with
            // foldTaxaBelow set, the latter is typically not even in the database's tree.
            TaxIdNode filedAt = foldUp(node, foldTaxaBelow, taxNodes);
            SmallTaxIdNode n = taxTree.getNodeByTaxId(filedAt.getTaxId());
            if (n == null) {
                return;
            }
            if (Rank.GENOME.ordinal() != n.getRankOrdinal()) {
                int pos = accessionEnd();
                int keyEnd = 1 + GenomeKeyTrie.genomeKeyLength(target, 1, pos);
                SmallTaxIdNode descendant = n.getDescendantWithName(target, 1, keyEnd);
                if (descendant != null) {
                    n = descendant;
                }
            }
            if (Rank.GENOME.ordinal() == n.getRankOrdinal()) {
                genomeNode = n;
                genomeCounts = map.get(n.getTaxId());
            }
        }

        /**
         * Returns the exclusive end index of the accession within {@link #target}, the same cut
         * {@code ReworkingStoreFastaReader} uses to derive a genome's key from its accession.
         *
         * @return the exclusive end index of the accession
         */
        private int accessionEnd() {
            int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
            if (pos < 0) {
                pos = size;
                while (pos > 0 && (target[pos - 1] == '\n' || target[pos - 1] == '\r')) {
                    pos--;
                }
            }
            return pos;
        }

        /**
         * Looks the k-mer up in the database and, if found, tallies it against the current contig's
         * genome at the position on its path the stored node sits at.
         *
         * @return whether the k-mer was tallied
         */
        @Override
        protected boolean handleStore(long kmer) {
            if (genomeNode == null || genomeCounts == null) {
                // No genome resolved for this contig - nothing to tally against.
                return false;
            }
            String storedTaxId = store.getLong(kmer, null);
            if (storedTaxId == null) {
                // K-mer is not in the database - it did not land anywhere.
                return false;
            }
            SmallTaxIdNode storedNode = taxTree.getNodeByTaxId(storedTaxId);
            if (storedNode == null) {
                return false;
            }
            int index = 0;
            SmallTaxIdNode n = genomeNode;
            while (n != null && n != storedNode) {
                n = n.getParent();
                index++;
            }
            if (n == null) {
                // storedNode is not on genomeNode's own path to the root - should not happen given the
                // lowest-common-ancestor update the database was built with, but guarded rather than
                // indexed out of bounds.
                if (getLogger().isWarnEnabled()) {
                    getLogger().warn("K-mer stored at " + storedTaxId + " is not on the path of genome "
                            + genomeNode.getTaxId() + ".");
                }
                return false;
            }
            // The array is the one pre-created in doMakeThis() and shared by every reader that meets this
            // genome - contigs of one genome may be read by several threads at once - so locking on it
            // yields per-genome mutual exclusion. The map itself is only read here, never modified.
            long[] counts = genomeCounts;
            synchronized (counts) {
                counts[index]++;
            }
            return true;
        }
    }
}
