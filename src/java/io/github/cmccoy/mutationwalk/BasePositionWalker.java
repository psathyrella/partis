package io.github.cmccoy.mutationwalk;

import java.io.FileWriter;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
//import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;

public abstract class BasePositionWalker extends LocusWalker {
    @Output
    public FileWriter out;

    //@ArgumentCollection
    //public StandardVariantContextInputArgumentCollection invrns = new StandardVariantContextInputArgumentCollection();
}
