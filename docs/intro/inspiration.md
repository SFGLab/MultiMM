# Inspiration

MultiMM started with a simple request. My friends in the lab are biologists. They had beautiful Hi-C maps, loop calls and compartment tracks, and they wanted one thing: *to see what the genome looks like in 3D*, quickly, without first having to learn polymer physics or read a pile of papers about force fields.

That request stayed with me. Over the years I have built several programs for 3D chromatin modelling. Each of them was, in the end, *a model*: a clever idea that works well for one question, at one scale, for the person who wrote it. MultiMM is meant to be something different. **It is a tool, not just a model.**

## What I wanted to build

My goal was a really **user-friendly** piece of software. In practice, that means a few things.

- **You should not need to know the physics to get a result.** Choose a level (a gene, a region, a chromosome or the whole genome), point MultiMM to your data, and run it. The defaults are chosen so that this works.
- **You should be able to model any chromatin scale.** The same program builds a 1000-bead picture of one gene and its surroundings, and a 200 000-bead picture of the entire nucleus. You do not switch software when your question grows.
- **Everything should be in one place.** All the forces (loops, compartments, Hi-C restraints, the nuclear envelope and more) and all their alternative forms are available as options. You can use a single one or combine them as you like.
- **It should be usable in different ways.** You can run it with one command, drive it from a small configuration file, or call it from your own Python code. It runs models, but it also checks them against your data and produces figures, so the whole path from data to picture is covered.

## Physics is a helper, not a barrier

Chromatin modelling is, in a good sense, *phenomenological*. We do not simulate every atom. We describe chromatin as a chain of beads and add simple rules that capture what the experiments show: beads cannot overlap, loops hold two points together, active and inactive regions like to meet their own kind. You do not need deep mathematics to use such rules, and many ideas in this documentation are as simple as "a spring" or "a magnet".

At the same time, I know that many people are scared of modelling. It looks like a locked room that only physicists can enter. I would like this documentation, and the software itself, to be a door that is open: if you are curious, you can go deeper and see exactly what each term does, with equations and pictures; if you only need a structure for your figure, you can stay on the surface and still get a sound result.

If you are a biologist, I hope MultiMM lets you test your ideas about genome folding in an afternoon. If you are a physicist or a computer scientist, I hope the [Model](../model/index.md) chapter gives you everything you need to understand, question and extend it.
