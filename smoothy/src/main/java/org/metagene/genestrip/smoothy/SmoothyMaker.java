package org.metagene.genestrip.smoothy;

import org.metagene.genestrip.GSMaker;

public class SmoothyMaker <P extends SmoothyProject> extends GSMaker<P> {
    /**
     * Creates a smoothy maker for the given project.
     *
     * @param project the project to create goals for
     */
    public SmoothyMaker(P project) {
        super(project);
    }
}
