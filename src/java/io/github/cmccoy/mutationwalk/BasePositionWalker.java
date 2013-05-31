package io.github.cmccoy.mutationwalk;

import java.io.PrintStream;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
//import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;

public abstract class BasePositionWalker extends LocusWalker {
    @Output
    public PrintStream out;

    //@ArgumentCollection
    //public StandardVariantContextInputArgumentCollection invrns = new StandardVariantContextInputArgumentCollection();
}
