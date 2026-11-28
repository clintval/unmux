wildcard_constraints:
    recipe="|".join(RECIPES),


for group in cookbook.fetch_groups():

    rule:
        name:
            cookbook.rule_name(group)
        output:
            [download(source.key, source.name) for source in group],
        log:
            logfile(cookbook.rule_name(group)),
        params:
            names=[(source.key, source.name) for source in group],
        run:
            cookbook.fetch_rule(CACHE, params, log)


rule tools:
    input:
        config["unmux"],
    output:
        TOOLS,
    log:
        logfile("tools"),
    params:
        lima=config["lima"],
    run:
        cookbook.tools_rule(CACHE, input, params, log)


rule versions:
    input:
        TOOLS,
    output:
        result("versions.json"),
    log:
        logfile("versions"),
    run:
        cookbook.versions_rule(CACHE, output, log)


rule check_recipe:
    input:
        run=work("{recipe}", "run.json"),
        script="scripts/checks/{recipe}.py",
    output:
        result("{recipe}.json"),
    log:
        logfile("{recipe}.check"),
    run:
        cookbook.check_rule(input, output, log)


rule render_recipe:
    input:
        results=result("{recipe}.json"),
        template="templates/{recipe}.md",
    output:
        result("{recipe}.md"),
    log:
        logfile("{recipe}.render"),
    run:
        cookbook.render_rule(input, output, log)


rule assemble:
    input:
        fragments=expand(result("{recipe}.md"), recipe=RECIPES),
        results=expand(result("{recipe}.json"), recipe=RECIPES),
        versions=result("versions.json"),
        order="config/recipes.yaml",
    output:
        result("recipes.md"),
    log:
        logfile("assemble"),
    run:
        cookbook.assemble_rule(input, output, log)


rule page:
    input:
        result("recipes.md"),
    output:
        "../recipes.md",
    log:
        logfile("page"),
    run:
        cookbook.page_rule(input, output, log)


rule check:
    input:
        result("recipes.md"),
    log:
        logfile("check"),
    run:
        cookbook.check_page_rule(input, log)
