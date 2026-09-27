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

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.smoothy.SmoothyProject;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Writes the per-genome path counts of {@link KMerGenomePathCountsGoal} to a CSV file, so that the
 * distributions p(n | m) can be checked and plotted.
 * <p>
 * The file has one row per genome and node on the genome's path to the root, in long format since the
 * paths differ in length. The columns are:
 * <ul>
 * <li>{@code genome taxid} and {@code genome name} - the genome's {@code GENOME} node, whose name is its
 * genome key (e.g. {@code NZ_ABCD} for a WGS assembly);</li>
 * <li>{@code position} - the index along the path, {@code 0} being the genome's own node;</li>
 * <li>{@code node taxid}, {@code node name} and {@code node rank} - the node at that position;</li>
 * <li>{@code kmers} - the genome's k-mer occurrences the database holds at that node;</li>
 * <li>{@code share} - {@code kmers} divided by the genome's total over its path, i.e. p(n | m). A
 * genome none of whose k-mers is in the database has no distribution, and its shares are left
 * empty.</li>
 * </ul>
 * Genomes are written in the order of the database's taxonomy tree, so the file is the same from one
 * run to the next.
 *
 * @param <P> the project type
 */
public class KMerGenomePathCountsCSVGoal<P extends SmoothyProject> extends FileGoal<P> {
    private static final DecimalFormat DF = new DecimalFormat("0.00000000", new DecimalFormatSymbols(Locale.US));

    private final ObjectGoal<Map<String, long[]>, P> pathCountsGoal;
    private final ObjectGoal<Database, P> dbGoal;

    /**
     * Creates the goal writing the path counts CSV file.
     *
     * @param project        the project this goal belongs to
     * @param key            the key identifying this goal
     * @param pathCountsGoal the goal providing the per-genome path counts
     * @param dbGoal         the goal supplying the loaded database, whose tree resolves the path nodes
     * @param deps           any further goals this goal depends on
     */
    @SafeVarargs
    public KMerGenomePathCountsCSVGoal(P project, GoalKey key, ObjectGoal<Map<String, long[]>, P> pathCountsGoal,
                                       ObjectGoal<Database, P> dbGoal, Goal<P>... deps) {
        super(project, key, Goal.append(deps, pathCountsGoal, dbGoal));
        this.pathCountsGoal = pathCountsGoal;
        this.dbGoal = dbGoal;
    }

    /**
     * The file sits in the project's {@code db} folder beside the database it was measured on, since
     * the counts are a property of that database and not of any sample. Its name is derived from the
     * project and goal names, as for every other output file.
     */
    @Override
    public List<File> getFiles() {
        return Collections.singletonList(getProject().getOutputFile(getProject().getDBDir(), getKey().getName(),
                null, null, GSProject.GSFileType.CSV, false));
    }

    @Override
    protected void makeFile(File file) throws IOException {
        Map<String, long[]> counts = pathCountsGoal.get();
        SmallTaxTree taxTree = dbGoal.get().getTaxTree();

        try (PrintStream ps = new PrintStream(file, StandardCharsets.UTF_8.name())) {
            ps.println("genome taxid;genome name;position;node taxid;node name;node rank;kmers;share;");
            for (SmallTaxIdNode genome : taxTree) {
                long[] genomeCounts = counts.get(genome.getTaxId());
                if (genomeCounts == null) {
                    continue;
                }
                long total = 0;
                for (long c : genomeCounts) {
                    total += c;
                }
                SmallTaxIdNode node = genome;
                for (int i = 0; i < genomeCounts.length && node != null; i++, node = node.getParent()) {
                    ps.print(genome.getTaxId());
                    ps.print(';');
                    ps.print(genome.getName());
                    ps.print(';');
                    ps.print(i);
                    ps.print(';');
                    ps.print(node.getTaxId());
                    ps.print(';');
                    ps.print(node.getName());
                    ps.print(';');
                    Rank rank = node.getRank();
                    ps.print(rank == null ? "" : rank.getName());
                    ps.print(';');
                    ps.print(genomeCounts[i]);
                    ps.print(';');
                    if (total > 0) {
                        ps.print(DF.format((double) genomeCounts[i] / total));
                    }
                    ps.println(';');
                }
            }
        }
    }
}
