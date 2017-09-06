#!/usr/bin/perl
use warnings;
use Data::Dumper;
use Getopt::Std;
use Config::Tiny;
use Time::HiRes;

#引数指定：初期個体数、各遺伝子頻度、交配様式の指定、時系列の出力ファイル
my %opts = ();
getopts("i:a:b:o:g:f:y:m:", \%opts);
my $StartIndividuals = $opts{i};
my $Alphas = $opts{a};
my $Betas = $opts{b};
my $Omegas = $opts{o};
my $Model = $opts{g} //= "GeneModel1" ;
my $MateModel = $opts{m} //= "MateModel1" ;
my $output = $opts{f};
my $year = $opts{y} //= "1";

#割合に変換
my $all = $Alphas + $Betas + $Omegas;
$Alphas = $Alphas / $all * 100;
$Betas = $Betas / $all * 100;
$Omegas = $Omegas / $all * 100;

#同一ディレクトリ内の設定ファイル読み込み
my $cabin = $0;
$cabin =~ s/omegaverse\.pl$//;
my $inicabin = $cabin . "omegaverse.ini";
my $init;
unless(-e $inicabin){
	print "[ERROR] Not Found '$inicabin' in '$cabin'\n";
	exit;
}else{
    $init = Config::Tiny -> read("$inicabin");
}

#各個体のプロファイル作成
my $Omegaverse;
#ID用の素数列読み込み
my $prime = $0;
$prime =~ s/omegaverse\.pl$//;
$prime = $prime . "primes.txt";
open(PR, "<$prime") or die;
my @primes = ();
my $c = 0;
while(<PR>){
    $c++;
    chomp;
    push(@primes, $_);
    last if $c == $StartIndividuals + 1;
}
#個体ごとの処理
for($i = 1; $i <= $StartIndividuals; $i++){
    my %Individual;
    $Individual{ID} = $primes[$i];
    $Individual{ParentID} = "0"; #初期親

    #ヒト性別遺伝子型
    my $XYtype = rand(100);
    if ($XYtype > 50){
        $Individual{XYtype} = "XX";
        $Individual{XYcode} = "1";
    }else{
        $Individual{XYtype} = "XY";
        $Individual{XYcode} = "0";
    };

    #オメガ遺伝子型
    my @ABOgene = ();
    for($j = 0; $j < 2; $j++){
        my $gene = rand(100);
        if($gene <= $Alphas){
            push(@ABOgene, "A"); #alpha
        }elsif($gene <= $Betas){
            push(@ABOgene, "B"); #beta
        }else{
            push(@ABOgene, "O"); #omega
        };
    }
    @ABOgene = sort{$a cmp $b} @ABOgene;
    $Individual{ABO} = join("", @ABOgene);

    #オメガ表現型(サブルーチン'Omegatype'で遺伝子モデルごとに決定)
    my %code = ("Alpha" => 2, "Beta" => 1, "Omega" => 0);
    my $omegatype = &Omegatype(\%Individual);
    $Individual{omegatype} = $omegatype;
    $Individual{omegacode} = $code{$omegatype};

    #年齢(世界統計2000年時点を採用、設定ファイルに記載)
    my $age = rand(100);
    if($age < $init->{population}->{under15}){ #15歳未満
        $age = int(rand(15)+0.5);
    }elsif($age >= (100-$init->{population}->{over65})){ #65歳以上
        my $lifespan = $init->{population}->{lifespan} - 65;
        $age = int(rand($lifespan)+0.5) + 65;
    }else{ #15~65歳
        my $lifespan = 50;
        $age = int(rand($lifespan)+0.5) + 15;
    }
    $Individual{age} = $age;

    #つがい(初期は無しとする？)
    my $matepair = 0;
    $Individual{matepair} = $matepair;
    $Individual{matenum} = $matepair;

    #状態(交配可能か、そうでないか)
    my $status = &Matestatus(\%Individual); #交配状態判定用ルーチンへ送る
    $Individual{status} = $status;

    #個体群に代入
    $Omegaverse->{$i} =  \%Individual;
}

#初期状態のカウント
print "---Inital State---\n";
&IndvidualStat();
print "------------------\n";
#print Dumper $Omegaverse;

#つがいの形成
#無作為に2個体抽出して番わせる
#人口千人比の婚姻率が5パーセント・出生率7パーセント前後になる程度に調整
#経過年数
for($k = 1; $k <= $year; $k++){
    #結婚処理
    $Omegaverse = &Mating($Omegaverse);
    #繁殖処理
    &Generation($Omegaverse);
    &IndvidualStat();
    #死亡・生存処理
}


#一年進ませて、平均余命から死亡判定

######################################

#個体群の統計値モニタリング
sub IndvidualStat {
    #カウント対象のkeyを作成(デフォルトはXYとABO)
    my $countkey = $_[0] //= "ABO|XYtype|omegatype";
    
    #各遺伝子型のカウント
    my %stats;
    foreach my $id (keys($Omegaverse)){
        foreach my $key (keys(%$Omegaverse{$id})){
           next unless ($key =~ /$countkey/);
           my $type = $$Omegaverse{$id}{$key};
           $stats{$type} //= 0;
           my $num = $stats{$type} + 1;
           $stats{$type} = $num;
        }
    }

    if($countkey =~ /ABO\|XYtype/){
        #モデルごとに表現型比、遺伝子型比を計算
        my $total = $stats{XX} + $stats{XY};
        print "Population size: $total\n";
        print "Sex:             ";
        #出力
        #sex
        foreach my $keys (sort keys(%stats)){
            next if ($keys =~ /A|B|O/);
            my $per = $stats{$keys} / $total * 100;
            $per = sprintf("%.2f", $per); #小数点下2ケタまで
            print "$keys: $per%\t";
        }
        print "\n";

        #omegatype
        print "Omegatype:       ";
        foreach my $keys (sort keys(%stats)){
            next if ($keys =~ /^..$/);
            my $per = $stats{$keys} / $total * 100;
            $per = sprintf("%.2f", $per); #小数点下2ケタまで
            print "$keys: $per%\t";
        }
        print "\n";

        #genotype
        print "Genotype:        ";
        foreach my $keys (sort keys(%stats)){
            next if ($keys =~ /XX|YY|XY|.{3,5}/);
            my $gene = $stats{$keys} / $total * 100;
            $gene = sprintf("%.2f", $gene);
            $stats{$keys} = $gene;
            print "$keys: $stats{$keys}% ";
            }
        print "\n";
    }
    #年齢分布
    my %agehist;
    if($countkey =~ /age/){
        my $total = 0;
    }
    #婚姻率分布
    if($countkey =~ /status/){
        print Dumper %stats;
    }
    if($countkey =~ /matepair/){
        print Dumper %stats;
    }

}
#ファイル出力
sub Outstat {
    open(OUT, ">$output");
        #header
        my @head = sort (keys(%$Omegaverse{$primes[0]}));
        print OUT "#";
        foreach my $head (@head){
            print OUT "$head\t";
        }
        print OUT "\n";
        foreach my $id (@primes){
            foreach my $elem (sort keys(%$Omegaverse{$id})){
                print OUT "$$Omegaverse{$id}{$elem}\t";
            }
            print OUT "\n";
        }
    close(OUT);
}
#オメガ表現型判定(設定モデルごと)
sub Omegatype {
    my $in = $_[0];
    my %input = %$in;
    my $type;
    #GeneModel1 (A>B>O, XY独立)
    if($Model eq "GeneModel1"){
        if($input{ABO} =~ /A/){
            $type = "Alpha";
        }elsif($input{ABO} =~ /B/){
            $type = "Beta";
        }else{
            $type = "Omega";
        }
    #GeneModel2 (A＝O > B, AO = B)
    }elsif($Model eq "GeneModel2"){
        if($input{ABO} =~ /A/){
            $type = "Alpha";
            $type = "Beta" if $input{ABO} eq "AO";
        }elsif($input{ABO} =~ /O/){
            $type = "Omega";
        }else{
            $type = "Beta";
        }
    }
    return($type);
}
#繁殖状態判定
sub Matestatus {
    my $in = $_[0];
    my %input = %$in;
    my $status = 1;
    #初潮前は繁殖不可(仮：男女共通にしておく)
    $status = 0 if ($input{age} < $init->{mating}->{menophania});
    #閉経後は繁殖不可
    $status = 0 if ($input{omegatype} =~ /O/ and $input{age} > $init->{mating}->{menopause});
    $status = 0 if ($input{XYtype} =~ /XX/ and $input{age} > $init->{mating}->{menopause});
    #つがいが確定している個体は繁殖不可
    $status = 0 unless ($input{matepair} eq 0);
    return($status);
}

#繁殖プロセス(繁殖可不可の判定→子の生成)
sub Mating {
    my $in = $_[0];
    my @pop = sort keys(%$in);
    my $population = @pop;
    my $mate_num = $population * $init->{mating}->{mating_per};
    my $mate_success = 0;
    my $birth_sucsess = 0;
    my $test_num = 0;
    #mate pair生成
    for($j = 1; $j <= int($mate_num); $j++){
        my $mate1 = int(rand($population)+0.5);
        my $mate2 = int(rand($population)+0.5);
        #自殖(同一ID)の場合は再代入
        while($mate1 == $mate2){
            $mate2 = int(rand($population)+0.5);
        }
        #$test_num++; #debug
        #statusに基づき交配可能判定
        next if $$in{$mate1}{status} eq 0;
        next if $$in{$mate2}{status} eq 0;

        #性別に基づいた子の生成可能判定
        my $materesult_omega = $$in{$mate1}{omegacode} + $$in{$mate2}{omegacode};
        my $materesult_XY = $$in{$mate1}{XYcode} + $$in{$mate2}{XYcode};
       
        #$test_num++; #debug
        #MateModel1(XX-XY, XX-XXかつ一方がA, XY-XYかつ一方がO, 同性でA-A, O-Oは除外)
        if($MateModel eq "MateModel1"){
            #同性かつアルファまたはオメガが無いケース,両方がアルファまたはオメガのケースを除外
            if($materesult_XY eq 0){ #XY-XYの場合
                #同性betaを除外
                next if($$in{$mate1}{omegacode} * $$in{$mate2}{omegacode} eq 1);
                #一方のみがOでないと子供は生まれない
                next if($materesult_omega =~ /0|3|4/);
            }elsif($materesult_XY eq 2){  #XX-XXの場合
                #同性betaを除外
                next if($$in{$mate1}{omegacode} * $$in{$mate2}{omegacode} eq 1);
                #一方がAでないと子供は生まれない
                next if ($materesult_omega =~ /0|1|4/);
            }
        #MateModel2
        }elsif($MateModel eq "MateModel2"){
            exit;
        #エラー処理
        }else{
            print "[ERROR]:'$MateModel' is not supported mating model\n";
            exit;
        }
        $mate_success++;

        #つがいの判定(A-O, およびB-B間ではつがい形成)
        if($$in{$mate1}{omegacode} + $$in{$mate1}{omegacode} eq 2){
            $$in{$mate1}{matepair} = $$in{$mate2}{ID};
            $$in{$mate1}{matenum} = $mate2;
            $$in{$mate2}{matepair} = $$in{$mate1}{ID};
            $$in{$mate2}{matenum} = $mate1;
        }
        #血縁度判定:素因数分解して共通する親の個数が2個
        #17/9/6 血縁度処理はいったん放置：
        #my $kindred = &Kindred($$in{mate1}{ID}, $$in{mate1}{ID}, $population); #サブルーチンへ投げる
        #next if $kindred < 3;
        
        #出生率調整
        my $birth_fate = rand(1);
        next if $birth_fate > $init->{mating}->{birth_rate};
        $birth_success++;
    }
    #結婚率出力
    my $mating_rate = $mate_success / $population * 1000;
    $mating_rate = sprintf("%.2f", $mating_rate);
    print "Mating rate: $mating_rate\n";
    return($in);

    #出生率出力
    #my $birth_rate = $birth_success / $population * 1000;
    #$birth_rate = sprintf("%.2f", $birth_rate);
    #print "Birth rate: $birth_rate\n";

    #print Dumper $mate_success / $test_num; 
}

#与えられた2数を素因数分解して、共通する因数の個数を返す
sub Kindred {
#test.plで開発中
}
#親の情報から子供の情報生成、omegaverseに追加
sub Generation {
    my $in = $_[0];
    my @pop = sort keys(%$in);
    my $population = @pop;
    foreach my $id (sort keys(%$in)){
        next if $$in{$id}{matepair} eq 0; #非婚は除外
        next if $$in{$id}{status} eq 2;   #出産済の場合は除外
        #親の情報取り込み
        my $parentid_1 = $$in{$id}{ID};
        my $pair_id = $$in{$id}{matenum};
        my $parentid_2 = $$in{$pair_id}{ID};
        #親の情報変更
        $$in{$id}{status} = 2;
        $$in{$pair_id}{status} = 2;

        #子供の情報生成 parentIDは血縁度の繁殖制限使わない限り無駄そう
        my %child;
        $child{ID} = $parentid_1 * $parentid_2;
        $child{ParentID} = $parentid_1 . "_" . $parentid_2;
        $child{age} = 0;
        $child{status} = 0;
        $child{matepair} = 0;
        $child{matenum} = 0;

        #子供の遺伝子型生成
        #XY
        my $XY_1 = $in{$id}{XYtype};
        my $XY_2 = $in{pair_id}{XYtype};
        my $gamete_xy_1 = substr($$in{$id}{XYtype}, int(rand(1)+1), 1);
        my $gamete_xy_2 = substr($$in{$pair_id}{XYtype}, int(rand(1)+1) ,1);
        my @gamete_xy = ($gamete_xy_1, $gamete_xy_2);
        @gamete_xy = sort{$a cmp $b} @gamete_xy;
        $child{XYtype} = join("",@gamete_xy);
        my %code_xy = ("XX" => 1, "XY" => 0, "YY" => -1);
        $child{XYcode} = $code_xy{$child{XYtype}};

        #ABO
        my $ABO_1 = $in{$id}{ABO};
        my $ABO_2 = $in{pair_id}{ABO};
        my $gamete_abo_1 = substr($$in{$id}{ABO}, int(rand(1)+1), 1);
        my $gamete_abo_2 = substr($$in{$pair_id}{ABO}, int(rand(1)+1) ,1);
        my @gamete_abo = ($gamete_abo_1, $gamete_abo_2);
        @gamete_abo = sort{$a cmp $b} @gamete_abo;
        $child{ABO} = join("", @gamete_abo);
        my $omegatype = &Omegatype(\%child);
        $child{omegatype} = $omegatype;
        my %code_abo = ("Alpha" => 2, "Beta" => 1, "Omega" => 0);
        $child{omegacode} = $code_abo{$omegatype};

        #個体群に代入
        $population++;
        $Omegaverse->{$population} = \%child;
    }
}
#各集団の年齢に1歳加える、交配ステータス処理(0→1、2→1)
#年齢および遺伝子型(YY)から一定確率で死亡判定、通し番号を整理して欠番をゼロに
sub Death {
#死亡者の相手のmate pairを空欄にする
#betaの場合はstatus解除(つがいでないため)
}
#####################################


#課題：配偶子生成と遺伝子型の優劣
#表現型ルール、独立・連鎖ルール・交配ルール
#個体間の血縁度の実装
#年齢・性別・オメガ性別
#つがい状態：交配可・交配不可
#年齢
#遺伝子型ルールと交配ルール
#初期集団は同じで、交配ルールだけ変えた試行を行う（forkするかしないかは別）
#交配ルーチン
#ターン制
#発生・死亡判定
#1.交配させるインスタンスを選択
#2.交配様式に従い交配の成否を決定（結果主義）
#3.交配集団に新規個体の追加
#####################################
#以下交配様式のサブルーチン
#交配可能年齢の処理
#
#XY性別、ABO性別に対しての処理
#配偶子の組み合わせ
#alpha-alpha/alpha-beta/alpha-omega/beta-beta/beta-omega/omega-omega/
#XY-XX/XY-XY/XX-XX
