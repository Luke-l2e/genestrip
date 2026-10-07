package org.metagene.genestrip.smoothy.goals;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.smoothy.SmoothyProject;

import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Writes the per-genome path counts of {@link KMerGenomePathCountsGoal} to a serialized Java object
 * file, so that later goals or analyses can reuse them without another pass over the RefSeq release.
 * <p>
 * The file holds the goal's result as is: a {@code HashMap<String, long[]>} keyed by a genome's tax id,
 * whose array holds the genome's k-mer occurrences per position on its path to the root, {@code 0}
 * being the genome's own node. The nodes themselves are not stored, since they follow from the tax
 * tree of the database the counts were measured on; {@link #load(File)} reads the map back.
 *
 * @param <P> the project type
 */
public class KMerGenomePathCountsSERGoal<P extends SmoothyProject> extends FileGoal<P> {
    private final ObjectGoal<Map<String, long[]>, P> pathCountsGoal;

    /**
     * Creates the goal writing the path counts SER file.
     *
     * @param project        the project this goal belongs to
     * @param key            the key identifying this goal
     * @param pathCountsGoal the goal providing the per-genome path counts
     * @param deps           any further goals this goal depends on
     */
    @SafeVarargs
    public KMerGenomePathCountsSERGoal(P project, GoalKey key, ObjectGoal<Map<String, long[]>, P> pathCountsGoal,
                                       Goal<P>... deps) {
        super(project, key, Goal.append(deps, pathCountsGoal));
        this.pathCountsGoal = pathCountsGoal;
    }

    /**
     * The file sits in the project's {@code db} folder beside the database it was measured on, as the
     * CSV file of {@link KMerGenomePathCountsCSVGoal} does. Its name is derived from the project and
     * goal names, as for every other output file.
     */
    @Override
    public List<File> getFiles() {
        return Collections.singletonList(getProject().getOutputFile(getProject().getDBDir(), getKey().getName(),
                null, null, GSProject.GSFileType.SER, false));
    }

    @Override
    protected void makeFile(File file) throws IOException {
        // Copied into a HashMap so that what is written is serializable whatever map the goal returns.
        HashMap<String, long[]> counts = new HashMap<>(pathCountsGoal.get());
        try (ObjectOutputStream oOut = new ObjectOutputStream(StreamProvider.getOutputStreamForFile(file))) {
            oOut.writeObject(counts);
        }
    }

    /**
     * Reads back path counts written by this goal.
     *
     * @param file the SER file to read
     * @return the per-genome path counts, keyed by the genome's tax id
     * @throws IOException            if the file cannot be read
     * @throws ClassNotFoundException if the serialized map cannot be resolved
     */
    @SuppressWarnings("unchecked")
    public static Map<String, long[]> load(File file) throws IOException, ClassNotFoundException {
        try (ObjectInputStream oIn = new ObjectInputStream(StreamProvider.getInputStreamForFile(file))) {
            return (Map<String, long[]>) oIn.readObject();
        }
    }
}
